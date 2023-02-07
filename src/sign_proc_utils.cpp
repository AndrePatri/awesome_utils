#include "include/awesome_utils/sign_proc_utils.hpp"

using std::tanh;

using namespace SignProcUtils;

//************* NumDiff *************//

NumDiff::NumDiff()
{

}

NumDiff::NumDiff(int n_jnts, double dt, int order)
    :_n_jnts{n_jnts}, _order{order}, _dt{dt}
{

    _k_n1 = Eigen::VectorXd::Zero(2); // 1st order
    _k_n2 = Eigen::VectorXd::Zero(3); // 2nd order
    _k_n3 = Eigen::VectorXd::Zero(4); // 3rd order
    _k_n4 = Eigen::VectorXd::Zero(5); // 4th order
    _k_n5 = Eigen::VectorXd::Zero(6); // 5th order
    _k_n6 = Eigen::VectorXd::Zero(7); // 6th order

    _k_n1 << 1.0,         -1.0; // 1st order
    _k_n2 << 3.0/2.0,     -2.0,   1.0/2.0; // 2nd order
    _k_n3 << 11.0/6.0,    -3.0,   3.0/2.0,   -1.0/3.0; //3rd order
    _k_n4 << 25.0/12.0,   -4.0,       3.0,   -4.0/3.0,    1.0/4.0; //4th order
    _k_n5 << 137.0/60.0,  -5.0,       5.0,  -10.0/3.0,  15.0/12.0,  -1.0/5.0; //5th order
    _k_n6 << 147.0/60.0,  -6.0,  15.0/2.0,  -20.0/3.0,   15.0/4.0,  -6.0/5.0,  1.0/6.0; //6th order

    _window_data = Eigen::MatrixXd::Zero(_n_jnts, _order + 1);

    if(_order == 1)
    {
        _k_n = _k_n1;
    }

    if(_order == 2)
    {
        _k_n = _k_n2;
    }

    if(_order == 3)
    {
        _k_n = _k_n3;
    }

    if(_order == 4)
    {
        _k_n = _k_n4;
    }

    if(_order == 5)
    {
        _k_n = _k_n5;
    }

    if(_order == 6)
    {
        _k_n = _k_n6;
    }

}

void NumDiff::add_sample(Eigen::VectorXd& sample)
{
  int sample_size = sample.size();

  if(sample_size != _n_jnts)
  {
      std::string exception = std::string("NumDiff::add_sample(): Trying to add a sample of size ") +
                              std::to_string(sample_size) + std::string(", which is different from ") +
                              std::to_string(_n_jnts) + std::string(", (number of joints) \n");

      throw std::invalid_argument(exception);
  }

  // shifting data to the right (discarting most remote sample, which is the
  // one on the extreme right)
  for (int i = _order - 1; i >= 0; i--)
  {
     _window_data.block(0, i + 1 , _n_jnts, 1) = _window_data.block(0, i, _n_jnts, 1);

  }

  _window_data.block(0, 0 , _n_jnts, 1) = sample; // assign most recent sample

}

void NumDiff::dot(Eigen::VectorXd& sample_dot)
{

  if(_order > 6)
  {
      std::string exception = std::string("NumDiff::dot(): Support only up to the 6-th order derivative estimation is provided.");

      throw std::invalid_argument(exception);
  }

  sample_dot = (_window_data * _k_n) / _dt; // estimating derivative with backward method
  // of order _order

}

//************* NumInt *************//

NumInt::NumInt()
{

}

NumInt::NumInt(int n_jnts, double dt, double T_horizon)
    :_n_jnts{n_jnts}, _dt{dt}, _T_horizon{T_horizon}
{
    _n_intervals = std::round(_T_horizon / _dt);

    _n_samples = _n_intervals + 1;

    _window_data = Eigen::MatrixXd::Zero(_n_jnts, _n_samples);

}

void NumInt::add_sample(Eigen::VectorXd& sample)
{
  int sample_size = sample.size();

  if(sample_size != _n_jnts)
  {
      std::string exception = std::string("NumInt::add_sample(): Trying to add a sample of size ") +
                              std::to_string(sample_size) + std::string(", which is different from ") +
                              std::to_string(_n_jnts) + std::string(", (number of joints) \n");

      throw std::invalid_argument(exception);
  }

  // shifting data to the right (discarting most remote sample, which is the
  // one on the extreme right)
  for (int i = _n_samples - 1; i > 0; i--)
  {
     _window_data.block(0, i, _n_jnts, 1) = _window_data.block(0, i - 1, _n_jnts, 1);

  }

  _window_data.block(0, 0 , _n_jnts, 1) = sample; // assign most recent sample

}

void NumInt::get(Eigen::VectorXd& sample_integral)
{
    sample_integral = Eigen::VectorXd::Zero(_n_jnts);

    for(int i = _n_intervals; i > 0; i--)
    { // we integrate all the data in the window
        sample_integral = sample_integral +
                ( _window_data.block(0, i, _n_jnts, 1) +
                  _window_data.block(0, i - 1, _n_jnts, 1) ) / 2.0 * _dt;
    }
}

//************* NumIntRt *************//

NumIntRt::NumIntRt()
{

}

NumIntRt::NumIntRt(int n_jnts, double dt)
    :_n_jnts{n_jnts}
{

    _num_int = NumInt(n_jnts, dt, dt); // we exploit the NumInt class over a window size of 2 samples

    _int_km1 = Eigen::VectorXd::Zero(_n_jnts);

    _int_k = Eigen::VectorXd::Zero(_n_jnts);

}

void NumIntRt::add_sample(Eigen::VectorXd& sample)
{

  _num_int.add_sample(sample);

  int sample_size = sample.size();

  if(sample_size != _n_jnts)
  {
      std::string exception = std::string("NumIntRt::add_sample(): Trying to add a sample of size ") +
                              std::to_string(sample_size) + std::string(", which is different from ") +
                              std::to_string(_n_jnts) + std::string(", (number of joints) \n");

      throw std::invalid_argument(exception);
  }

}

void NumIntRt::get(Eigen::VectorXd& sample_integral)
{

    _num_int.get(_int_k);

    _int_k = _int_km1 + _int_k;

    sample_integral = _int_k;

    _int_km1 = _int_k;

}

//************* MovAvrgFilt *************//

MovAvrgFilt::MovAvrgFilt()
{

}

MovAvrgFilt::MovAvrgFilt(int n_jnts, double dt, int window_size)
    :_n_jnts{n_jnts}, _window_size{window_size}, _samples_dt{dt}
{

    if (_window_size <= 1)
    {
        _window_size = 2;
    }

    _cutoff_freq = _magic_number/_samples_dt * 1.0 / std::sqrt(std::pow(_window_size, 2) - 1);

    _window_data = Eigen::MatrixXd::Zero(_n_jnts, _window_size);

}

MovAvrgFilt::MovAvrgFilt(int n_jnts, double dt, double cutoff_freq)
    :_n_jnts{n_jnts}, _cutoff_freq{cutoff_freq}, _samples_dt{dt}
{

    double nominal_window_size = std::sqrt(std::pow(_magic_number / (_samples_dt * _cutoff_freq), 2) + 1);
    _window_size = std::round(nominal_window_size);
    if (_window_size <= 1)
    {
        _window_size = 2;
    }
    _window_data = Eigen::MatrixXd::Zero(_n_jnts, _window_size);

}

void MovAvrgFilt::add_sample(Eigen::VectorXd& sample)
{

  int sample_size = sample.size();

  if(sample_size != _n_jnts)
  {
      std::string exception = std::string("MovAvrgFilt::add_sample(): Trying to add a sample of size ") +
                              std::to_string(sample_size) + std::string(", which is different from ") +
                              std::to_string(_n_jnts) + std::string(", (number of joints) \n");

      throw std::invalid_argument(exception);
  }

  // shifting data to the right (discarting most remote sample, which is the
  // one on the extreme right)
  for (int i = _window_size - 2; i >= 0; i--)
  {
     _window_data.block(0, i + 1 , _n_jnts, 1) = _window_data.block(0, i, _n_jnts, 1);

  }

  _window_data.block(0, 0 , _n_jnts, 1) = sample; // assign most recent sample

  if (_is_first_run)
  { // if it's the first time we add a sample
    // we fill the empty part of the window with
    // replicae of the added sample.
    // after _window_size samples, the filter is at
    // regime.

    for (int i = 1; i < _window_data.cols(); i++)
    {
        _window_data.block(0, i, _n_jnts, 1) = sample;
    }

    _is_first_run = false;

  }

}

void MovAvrgFilt::get(Eigen::VectorXd& filt_sample)
{

    filt_sample = 1.0 / _window_size * _window_data.rowwise().sum();

}

void MovAvrgFilt::get_cutoff_freq(double& cutoff_f)
{
    cutoff_f = _cutoff_freq;
}

void MovAvrgFilt::get_window_size(int& window_size)
{
    window_size = _window_size;
}

//************* SignWithMem *************//

SignWithMem::SignWithMem()
{

}

SignWithMem::SignWithMem(double signal_3sigma,
                         double tanh_coeff)
    :_signal_3sigma{signal_3sigma}, _tanh_coeff{tanh_coeff}
{

    _tanh_thresh = tanh(_tanh_coeff * signal_3sigma);

}

double SignWithMem::approx_sign(double value)
{
    double approx_sign = tanh(_tanh_coeff * value);

    return approx_sign;
}

void SignWithMem::sign_with_memory()
{

    double sign_approximation = approx_sign(_value);

    if (sign_approximation > _tanh_thresh)
    {
        _sign = 1;
    }
    if (sign_approximation <= _tanh_thresh && sign_approximation >= -_tanh_thresh)
    { // uncertainty region --> retain previous value

        _sign = _previous_sign;
    }
    if (sign_approximation < -_tanh_thresh)
    {
        _sign = -1;
    }

    _previous_sign = _sign;
}

int SignWithMem::sign(double value)
{
    _value = value; // assign value

    sign_with_memory(); // compute sign

    return _sign;

}

//************* SmoothTanhSign *************//

SmoooothSign::SmoooothSign()
{

}

SmoooothSign::SmoooothSign(double signal_3sigma,
                           int alpha,
                           double beta,
                           bool use_threshold)
    :_signal_3sigma{signal_3sigma}, _alpha{alpha}, _beta{beta}
{

    if(alpha < 1)
    {
        std::string exception = std::string("SmoooothSign::SmoooothSign(): alpha should more or equal to 1");

        throw std::invalid_argument(exception);
    }

    _k = atanh(_beta) / (abs(_signal_3sigma) * _alpha);

    _threshold = _signal_3sigma;

    _use_threshold = use_threshold;

}

SmoooothSign::SmoooothSign(double signal_3sigma,
                           int alpha,
                           bool use_threshold)
    :_signal_3sigma{signal_3sigma}, _alpha{alpha}
{

    if(alpha < 1)
    {
        std::string exception = std::string("SmoooothSign::SmoooothSign(): alpha should more or equal to 1");

        throw std::invalid_argument(exception);
    }

    _k = atanh(_beta) / (abs(_signal_3sigma) * _alpha);

    _threshold = _signal_3sigma;

    _use_threshold = use_threshold;

}

double SmoooothSign::smooooth_sign(double value)
{
    if (_use_threshold && abs(value) <= _threshold)

        return 0.0;

    else
    {
        return tanh(_k * value);
    }

}

double SmoooothSign::sign(double value)
{

    return SmoooothSign::smooooth_sign(value);

}


























