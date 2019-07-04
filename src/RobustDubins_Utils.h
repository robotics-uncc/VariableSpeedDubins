#ifndef ROBUST_DUBINS_UTILS_H
#define ROBUST_DUBINS_UTILS_H

#include<vector>
#include<string>
#include<assert.h>

typedef std::vector<double> vd;
typedef std::vector<int> vi;

namespace RobustDubins {


void DubinsChain( double turnRadius, double maxSpacing, vd xPts, vd yPts, vd hRadPts, 
                  vd & xPath, vd & yPath, vd & hRadPath );

// determine the pathClass (BSB or BBB) from the pathType 
// (LSL,LSR,RSL,RSR,LRL,RLR)
std::string DubinsPathClass(const std::string & pathType);

// determine the orientation of each segment given a pathType
// e.g., LSR returns: 
// ki = 1 (Left init. turn)
// km = 1 (Straight mid. seg)
// kf = -1 (Right final seg.)
bool DubinsPathCurvatureSigns(const std::string & pathType, 
                              double & ki, double & km, 
                              double & kf);

// Note: 
// aUnsigned is arc-length of the first arc 
// bUnsigned can be (BBB) the arc-length of the second arc 
//           or     (BSB) the arc-length of the middle straight segment
// cUnsigned can be (BBB) the arc-length of the third arc 
//           or     (BSB) the arc-length of the second arc 
// the endpoint vector returned must be initialized with size 3 

// returns an intermediate point along the path, given the arclength and the 
// three a,b,c params defining the path 
bool computeDubinsPoint(const std::string & pathType, 
                           const double & aUnsigned, 
                           const double & bUnsigned, 
                           const double & cUnsigned,
                           const double & xInit,
                           const double & yInit,
                           const double & hInit,
                           const double & radius,
                           const double & arclength,
                           vd & endpoint);

// returns endpoint of either a BSB or BBB path 
bool computeDubinsEndpoint(const std::string & pathType,
                           const double & aUnsigned, 
                           const double & bUnsigned, 
                           const double & cUnsigned,
                           const double & xInit,
                           const double & yInit,
                           const double & hInit,
                           const double & R,
                           vd & endpoint);


// returns endpoint of a BBB path 
bool computeDubinsBBBendpoint(const std::string & pathType, 
                              const double & aUnsigned, 
                              const double & bUnsigned, 
                              const double & cUnsigned,
                              const double & xInit,
                              const double & yInit,
                              const double & hInit,
                              const double & R,
                              vd & endpoint);

// returns endpoint of a BSB path 
bool computeDubinsBSBendpoint(const std::string & pathType, 
                              const double & aUnsigned, 
                              const double & bUnsigned, 
                              const double & cUnsigned,
                              const double & xInit,
                              const double & yInit,
                              const double & hInit,
                              const double & R,
                              vd & endpoint);

} // namespace

#endif
