#include <cmath>
#include "utils/test/utils.h"
#include "math/include/matrix.h"
#include "math/include/constants.h"
#include "utils/include/debug.h"

namespace feasst {

TEST(Matrix, axis_angle) {
  RotationMatrix mat;
  Position axis;
  axis.set_vector({2., 0., 0.});
  mat.axis_angle(axis, 90);
  EXPECT_NEAR(mat.determinant(), 1., NEAR_ZERO);
  EXPECT_EQ(mat.multiply(axis).coord(), axis.coord());
  Position point1;
  point1.set_vector({0., 2., 0.});
  Position point2;
  point2.set_vector({0., 0., 2.});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));
  mat.axis_angle(axis, -135);
  point2.set_vector({0., -std::sqrt(2.), -std::sqrt(2.)});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));

  // For a rotation matrix, the inverse should be equal to the transpose
  RotationMatrix mat_t = mat, mat_inv = mat;
  mat_t.transpose();
  mat_inv.invert();
  EXPECT_TRUE(mat_t.is_equal(mat_inv));
  EXPECT_FALSE(mat.is_equal(mat_inv));
}

TEST(Matrix, 2d) {
  RotationMatrix mat;
  Position axis;
  axis.set_vector({0., 0.});
  mat.axis_angle(axis, 90);
  EXPECT_NEAR(mat.determinant(), 1., NEAR_ZERO);
  Position point1;
  point1.set_vector({2., 0.});
  Position point2;
  point2.set_vector({0., 2.});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));
  mat.axis_angle(axis, -135);
  point2.set_vector({-std::sqrt(2.), -std::sqrt(2.)});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));

  // For a rotation matrix, the inverse should be equal to the transpose
  RotationMatrix mat_t = mat, mat_inv = mat;
  mat_t.transpose();
  mat_inv.invert();
  EXPECT_TRUE(mat_t.is_equal(mat_inv));
  EXPECT_FALSE(mat.is_equal(mat_inv));
}

TEST(Matrix, multiply) {
  Matrix A({{1, 2, 3}, {4, 5, 6}});
  Matrix B({{7, 8}, {9, 10}, {11, 12}});
  Matrix C = A.multiply(B);
  EXPECT_EQ(C.value(0, 0), 58);
  EXPECT_EQ(C.value(0, 1), 64);
  EXPECT_EQ(C.value(1, 0), 139);
  EXPECT_EQ(C.value(1, 1), 154);
}

TEST(Matrix, is_identity) {
  Matrix mat({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
  EXPECT_TRUE(mat.is_identity());
}
TEST(RotationMatrix, euler_x) {
  RotationMatrix mat;
  mat.set_size(3, 3);
  Position x({0.6, 0.6, 0.6});
  mat.euler_x(PI/2., 0., 0.);
  Position x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(1), -0.6, NEAR_ZERO);
  mat.euler_x(0, PI/2., 0.);
  x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(2), -0.6, NEAR_ZERO);
  mat.euler_x(0., 0., PI/2.);
  x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(1), -0.6, NEAR_ZERO);
  mat.euler_x(PI, 0., 0.);
  x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(0), -0.6, NEAR_ZERO);
  
  // check inverse rotation
  Matrix transpose = mat;
  transpose.transpose();
  Matrix identity = mat.multiply(transpose);
  EXPECT_TRUE(identity.is_identity());
}

}  // namespace feasst
