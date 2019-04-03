// Creates quadrature rule on an edge
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef QUAD1D_HPP
#define QUAD1D_HPP
#include<stddef.h>
namespace HArDCore2D {  // forward dec
class LegendreGauss;
}

// \addtogroup Quadratures
//@{


namespace HArDCore2D {

class QuadRuleEdge {
    static constexpr size_t max_doe = 20;

public:
    QuadRuleEdge(size_t doe, bool warn);

    size_t nq();
    double xq(size_t i);
    double yq(size_t i);
    double wq(size_t i);
    void setup(double xV[], double yV[]);

private:
    size_t _npts;
    LegendreGauss* _rule;
    double* _w;
    double* _xr;
    double* _yr;
	double _length;
};
class LegendreGauss {
public:
    LegendreGauss(size_t doe);  //, double x[] double y[], double w[]);
    void sub_rule_01();
    void sub_rule_02();
    void sub_rule_03();
    void sub_rule_04();
    void sub_rule_05();
    void sub_rule_06();
    void sub_rule_07();
    void sub_rule_08();
    void sub_rule_09();
    void sub_rule_10();
	void sub_rule_11();
    size_t npts();
    double wq(size_t i);
    double tq(size_t i);

private:
    size_t _doe;
    size_t _npts;
    double* _t;
    double* _w;
};
//@}
}

#endif /* QUAD1D_HPP */
