#include <Eigen/Core>
#include <Eigen/Geometry>

#include <string>

namespace Eigen
{
    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    typedef Eigen::Matrix<double, 3, 1> Vector3d;
}

namespace ContactEstUtils
{

  class MomentumBasedFObs
  {
      public:

          MomentumBasedFObs();

      private:

          std::string _urdf_path;

  };

}
