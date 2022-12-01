#include <Eigen/Core>
#include <Eigen/Geometry>

#include <string>

#include "awesome_utils/sign_proc_utils.hpp"
#include "awesome_utils/model_interface.hpp"

using namespace SignProcUtils;
using namespace ModelInterface;

namespace Eigen
{
    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    typedef Eigen::Matrix<double, 3, 1> Vector3d;
}

namespace ContactEstUtils
{
  typedef Eigen::MatrixXd MatrixXd;
  typedef Eigen::VectorXd VectorXd;

  class MomentumBasedFObs
  {
      public:

          MomentumBasedFObs();

          MomentumBasedFObs(Model::Ptr model_ptr, double data_dt);
          MomentumBasedFObs(Model::Ptr model_ptr, double data_dt, double bandwidth);

          void update(std::string contact_framename);

          void get_tau_obs();
          void get_f_est(); // get force estimate
          void get_w_est();

      private:

          Model::Ptr _model_ptr;

          int _nv = -1.0;

          double _dt = -1.0;

          double _BW_red_factor = 0.707; // attenuation of the signal
          // of 3dB (definition of bandwidth)

          double _bandwidth = 0.0;
          double _k = 0.0;

          std::string _current_cont_frame;

          MatrixXd _K;

          NumInt _integrator_tau;
          NumInt _integrator_g;
          NumInt _integrator_C_T_q_dot;

          VectorXd _tau_c_obs; // observed joint contact efforts

  };

}
