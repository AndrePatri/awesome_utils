#include <Eigen/Core>
#include <Eigen/Geometry>

#include <string>

#include "include/awesome_utils/sign_proc_utils.hpp"
#include "include/awesome_utils/model_interface.hpp"

using namespace SignProcUtils;
using namespace ModelInterface;

//namespace Eigen
//{
//    typedef Eigen::Matrix<double, 6, 1> Vector6d;
//    typedef Eigen::Matrix<double, 3, 1> Vector3d;
//}

namespace ContactEstUtils
{
    typedef Eigen::MatrixXd MatrixXd;
    typedef Eigen::VectorXd VectorXd;

    /// \brief Class to perform a momentum-based force estimation.
    ///
    /// Consider the following expression of the rigid-body dynamics:
    ///
    /// (1.1) B(q) * a + C(q, v) * q_dot + g = tau + tau_c
    ///
    /// where
    /// - B is the joint-space inertia matrix
    /// - C is the matrix of Coriolis and centrifugal terms (recall that (2) B_dot = C + C^T)
    /// - g is the joint-space gravitational vector
    /// - tau is a vector of joint-space effort measurements
    /// - tau_c is a "disturbance" vector, which can be used
    ///   to estimate the contact forces, reflected to the joints.
    ///
    /// The dynamics of the generalized momentum is given by
    ///
    /// (1.2) p_dot = C^T*v - g + tau + tau_c
    ///
    /// The dynamics of the (linear) observer is given by
    ///
    /// (3) y_dot = K * (p_dot - tau + g - C^T * v - y)
    ///
    /// which is an asymptotically stable dynamics when is
    /// K positive definite square matrix (usually is chosen as a diagonal matrix).
    ///
    /// - p is the joint-space momentum, given by B(q) * v
    ///
    /// Equation (3) can be obtained easily by considering that
    ///
    /// (4) p_dot = B_dot * v + B * a
    ///
    /// substituting (1) and (2) into (3) and extracting tau_c as the desired variable
    /// to be observed, one obtains the observer dynamics (3).
    ///
    /// How to use (3)? It can be easily integrated (numerically) to obtain an estimate of tau_c
    /// with some noise rejection properties: the higher the gains, the lower the rejection
    /// of noise, the faster the convergence.
    /// In particular, supposing a matrix of diagonals values, the bandwidth of this observer
    /// can be approximated as (5) BW = - k/log(1 - 0.707)
    ///
    /// Clearly, the quality of the estimate depends on the quality of the measurements (accuracy, noise, etc..)
    /// but also on the accuracy up to which the inertial properties of the system are known
    ///
    /// a numerical implementation to obtain y at the current time (y_k) can be implemented integrating both sides of (3)
    /// over a sample interval:
    ///
    /// int_0^{h} [ y_dot * d_t] = K * int_0^{h} [ (p_dot - tau + g - C^T * v - y) * dt] --->
    ///
    /// y_k - y_km1 = K * (p_k - p_km1 + int_0^{h} [ g - tau - C^T * v - y ] * dt)
    ///
    /// let us approximate  int_0^{h} [ y * dt] as (y_k + y_km1)/2.0 * h (trapezoidal integration).
    ///
    /// Rearranging,  we obtain
    ///
    /// (6) (I + dt/2.0 * K) * y_k = (I - dt/2.0 * K) * y_km1 + K * (p_k - p_km1 + int_0^{dt} [ g - tau - C^T * v - y ] * dt)
    ///
    /// The term " int_0^{dt} [ g - tau - C^T * v - y ] * dt " can be simply approximated using again trapezoidal integration.
    ///
    /// Inverting (6) w.r.t. y_k gives the update equation for y, i.e. the observer of the residual joint torques
    ///
    /// Note that (6) does not require the differentiation of q_dot and is hence less prone to noise than other possible model-based
    /// observer implementations.
    ///
    /// To obtain the the contact wrenches consider that:
    /// y \approx tau_c = sum_0^{n_c - 1} J^T_{i} * w_i
    /// where n_c is the number of contacts, w_i is the wrench on the i-th contact and J_i is the analytical
    /// jacobian of the i-th contact (each contact is described by a frame).
    /// we can concatenate everything into
    ///
    /// y \approc tau_c = J^T_{tot} * W
    /// where
    /// J_{tot} = [J_1; J_2; ....; J_k; ....; J_{n_c - 1}]
    /// and W = [w_1; w_2; ...; w_k; ....; w_{n_c - 1}]
    ///
    /// Given a regularization W_{lambda},
    /// and a regularization diagonal matrix Lambda,
    /// this problem can be approached in a least-square error sense
    /// by solving the system
    /// [J^T_{tot}; sqrt(Lambda)] * W = [tau_c; sqrt(Lambda) * W_lambda]
    /// which can be synthetically written as
    /// A * W = b
    /// (recall that the least square regression problem is TBsolved is
    /// min_{W_c}{ 1/2 * ||A * W_c - b||^2 }
    /// )
    ///
    /// Side note on the solution to y_dot = K (tau_c - y) in scalar case
    /// with tau_c = A * sin(omega * t) and omega = 2 * PI * f_tau_c:
    /// y(t) = (y0 + A * k * w^2/ (k^2 + w^2)) * e^{-k*t} +
    ///        + A * k / sqrt(k^2 + w^2) * sin(omega * t + phi), with phi = atan(omega/k)
    ///
    /// which means the bandwidth of the observer can be computed as
    /// 2 * pi * f_bw = k (recall that the bandwidth is defined as a reduction of the amplitude of
    /// the frequency response of -3dB = 20 * log(1/sqrt(2)))

    class MomentumBasedFObs
    {
    public:

        typedef std::weak_ptr<MomentumBasedFObs> WeakPtr;
        typedef std::shared_ptr<MomentumBasedFObs> Ptr;
        typedef std::unique_ptr<MomentumBasedFObs> UniquePtr;

        typedef Matrix<double, 6, 1> Reg6D;

        MomentumBasedFObs();

        MomentumBasedFObs(Model::Ptr model_ptr, double data_dt,
                          std::vector<std::string> contact_framenames,
                          double bandwidth = 10.0,
                          double lambda = 1.0, bool regularize_f = true,
                          std::vector<int> selector = std::vector<int>{0, 1, 2, 3, 4, 5});
        MomentumBasedFObs(Model::Ptr model_ptr, double data_dt,
                          std::vector<std::string> contact_framenames,
                          double bandwidth = 10.0,
                          Reg6D lambda = Eigen::MatrixXd::Identity(6, 1), bool regularize_f = true,
                          std::vector<int> selector = std::vector<int>{0, 1, 2, 3, 4, 5});

        void update();

        void set_active_contacts(std::vector<int> contact_indeces); // activates the contacts given by
        // contact_indeces

        void get_tau_obs(VectorXd& tau_c); // get contact joint efforts estimate
        void get_w(Eigen::VectorXd& w);
        void get_w_est_at(std::string contact_framename, Model::Wrench& w_c); // get the wrench estimate at a specific contact frame
        void get_w_est_at(int contact_index, Model::Wrench& w_c);
        void get_f_est_at(std::string contact_framename, Model::Force3D& f_c); // get force estimate
        void get_f_est_at(int contact_index, Model::Force3D& f_c); // get force estimate
        void get_t_est_at(std::string contact_framename, Model::Torque3D& t_c); // get wrench estimate
        void get_t_est_at(int contact_index, Model::Torque3D& t_c);
        void get_contact_framenames(std::vector<std::string> names);
        void get_contact_indeces(std::vector<int> indeces);

    private:

        Model::Ptr _model_ptr;

        int _nv = -1, // dimension of the dynamics of the system (equal to the size of the generalized velocity v)
            _nc = 0;

        bool _regularize_f = false;// will use previous solution of f_c to regularize the new solution
        // instead of using always a constant value

        double _dt = -1.0;

        Reg6D _lambda; // regularization vector
        // for the contact wrench of each contact

        double _BW_red_factor = 0.70711; // attenuation of -3dB ( = 20 * log10 (1/sqrt(2)) )
        double _PI =  2 * std::acos(0);
        double _bandwidth = 10.0; // [Hz]
        double _k = 12.3; // [1/s] (or [Hz])

        std::map<std::string, int> _contact_map;

        std::vector<int> _selector{0, 1, 2, 3, 4, 5}; // used to select partial estimation of components
        std::vector<int> _deselector{0, 1, 2, 3, 4, 5};

        std::vector<std::string> _contact_framenames,
                                 _valid_contact_names; // results of filtering _contact_framenames

        std::vector<int> _contact_indeces; // this holds the order in which contacts where
        // added (by default in the same order of contact_framenames, so for now not so useful)

        std::vector<bool> _active_contacts;

        MatrixXd _K;
        MatrixXd _Skp1, _Skp1_inv, _Sk; // state transition matrices for the discrete integration
        // of the observer dynamics

        NumInt _integrator;

        VectorXd _tau_c_k; // observed joint contact efforts

        VectorXd _p_km1; // last joint-space momentum at k - 1 instant (i.e. previous)

        VectorXd _v, _tau, _g, _p,
                 _to_be_integrated, _integral;
        MatrixXd _C;

        MatrixXd _A,
                 _A_lambda;
        VectorXd _b, _b_lambda;

        Eigen::VectorXd _W; // estimated contact wrenches ([6 x 1 -> linear + angular] * 6)
        Eigen::VectorXd _W_reg; // regularization vector for W
        Model::Wrench _w_buff; // used only as a temporary holder

        Model::SpatialJac _J_buffer; // used only as a temporary holder
        Eigen::MatrixXd _J_c_tot; // vertical contactenation
        // of all (generalized) contact jacobians

        void assign_regression_matrices(); // assign values to the regression matrix A and vector b
        void get_contact_jacobians(); // gets and assign jacobian from all active contact frames
        void process_contactnames();// we process input contact names and eliminate non valid ones (
        // i.e. frames that do not exist in the model)
        void process_selector(); // processes the selector given to the constructor in a suitable way
        void compute_bandwidth();

        void compute_tau_c(); // computes the observed value of tau_c, i.e. the residual joint efforts

        void apply_component_selector(Model::Wrench& vector);
        void apply_component_selector(Model::SpatialJac& J);

    };

}
