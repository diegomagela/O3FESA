#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <vector>

struct point_weight
{
    const std::vector<double> points;
    const std::vector<double> weights;
};

static const std::vector<double> x1{0.0};
static const std::vector<double> w1{2.0};

static const std::vector<double> x2{0.57735, -0.57735};
static const std::vector<double> w2{1.0, 1.0};

static const std::vector<double> x3{0.0, 0.77459, -0.77459};
static const std::vector<double> w3{0.88889, 0.55556, 0.55556};

static const std::vector<double> x4{0.33998, -0.33998, 0.86113, -0.86113};
static const std::vector<double> w4{0.65214, 0.65214, 0.34785, 0.34785};

static const std::vector<point_weight> point_weight_vec = {{x1, w1},
                                                           {x2, w2},
                                                           {x3, w3},
                                                           {x4, w4}};

class Quadature
{
public:
    Quadature(const std::size_t n_points) : points_(point_weight_vec.at(n_points - 1).points),
                                            weights_(point_weight_vec.at(n_points - 1).weights){};

    // Selectors

    inline std::vector<double> get_points() const { return points_; }
    inline std::vector<double> get_weights() const { return points_; }

private:
    std::vector<double> points_;
    std::vector<double> weights_;
};

#endif // QUADRATURE_H
