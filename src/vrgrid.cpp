#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "csv.h"
#include <cstddef>
#include <cstring>
#include <random>
#include <cmath>
#include <vectorgraph.h>
#include <barycenterlayout.h>
#include <voronoy.h>
#include <wulipvector.h>
#include <incrementaldelaunay.h>
#include <polygon.h>
#include <Rectangle.h>
#include <map>
#include <unordered_set>

#include <numeric>
#include <chrono>
#include <thread>
// #include <graphvizfunc.h>
// #include <utils.h>
// #include <pythonLauncher.h>
// gui includes
#define GL_GLEXT_PROTOTYPES
#define EGL_EGLEXT_PROTOTYPES
extern "C" {
#include <GL/gl.h>
#include <GL/glut.h>
}
#include <Graph2DRenderer.h>
#include <Animation.h>
#include <fatum.h>
// /////////////////////
// Rendering objects
static wlp::Fatum* dbVis = 0;
static int WIDTH = 800;
static int HEIGHT = 800;
static int lastX=-1, lastY=-1;
static int mouseState;

bool useIterativeBorder = false; // False for spaceFilling before snap
float sensibility = 1e-2;
// float sensibility = 5e-3;
float gapSensibility = 5e-2;
bool useFixedIterationNumber = false;
// double maxRange = .3;
double maxRange = 100;
double maxRangeSnap = 2;
int snapCheckFreq = 8;
double forceSnapRatio = 0;
// double forceSnapRatio = 0.95;
int iterationInverseFactor = 2;
std::string outputFilename;

using namespace std;
using namespace wlp;
// Struct used in rendering
std::vector<Vec2d> nodePosOrig;

// classId == -1 => is border node
// classId == -2 => is dummy
// classId == -3 => is space filler

std::chrono::steady_clock::time_point beginC;

double lerp(double a, double b, double f) {
    return a + (b - a) * f;
}

wlp::Color getColorFromPosition(double x, double y, int targetResolution) {
    double xFactor = x/targetResolution + 0.5;
    double yFactor = y/targetResolution + 0.5;

    double r = lerp(0., 1., 1-xFactor);
    // double r = 1-xFactor*2;
    double g = lerp(0., 1., yFactor);
    double b = lerp(0., 1., xFactor);
    // double b = xFactor * 2;
    // double g = yFactor;
    return wlp::Color(int(r*255), int(g*255), int(b*255));

}

vector<Color> jetColormap{
    Color(0, 0, 131),
    Color(0, 60, 170),
    Color(5, 255, 255),
    Color(255, 255, 0),
    Color(250, 0, 0),
    Color(128, 0, 0)
};

Color lerpColor(Color c1, Color c2, float l) {
    if (l > 1) {
        l = 1;
    }

    if (l < 0) {
        l = 0;
    }
    int varR = c2.r() - c1.r();
    int varG = c2.g() - c1.g();
    int varB = c2.b() - c1.b();
    int vR = l * varR;
    int vG = l * varG;
    int vB = l * varB;
    return Color(c1.r() + vR, c1.g() + vG, c1.b() + vB);
}

Color getColorFromValue(double v, vector<Color> &colormap) {
    if (v > 1 || v < 0) {
        return Color(0, 0, 0);
    }

    v = max(v, 0.);
    double step = 1. / (colormap.size() - 1);
    int istep = 0;
    while (v > (istep * step)) {
        ++istep;
    }

    --istep;

    double vmin = istep * step;
    double vmax = (istep + 1) * step;
    Color c1 = colormap[istep];
    Color c2 = colormap[istep + 1];
    vmax = vmax - vmin;
    v = v - vmin;
    v = v / vmax;
    return lerpColor(c1, c2, v);
}

wlp::Color getColorFromClassId(int classId) {
    float v = (classId * 1.0 - 5.0) / 9.0 ;
    return getColorFromValue(v, jetColormap);
}

/**
 * @brief The SampleMark struct allows displaying "pxels" in the 2d space using position and class of samples
 */
struct SampleMark {
    SampleMark(wlp::Fatum *dbvis, const Vec2d &position, int classId, const int targetResolution):dbvis(dbvis) {
        dbvis->defaultMark().color(Color::Black)
                .size(Size(1, 1, 1))
                .shape(wlp::Shape::TRIANGLE);
        wlp::Shape s = wlp::Shape::SQUARE;
        wlp::Mark &mref = dbvis->addMark()
                .position(position.x(), position.y(), 0)
                .size(Size(.5,.5,0.))
                .shape(s);
        if (classId == -1) {
            mref.size(Size(0.5,.5,0.)).color(Color::Red)
            .shape(wlp::Shape::CIRCLE);            
        } else if (classId == -2) {
            mref.size(Size(0.5,.5,0.)).color(Color::White)
            .shape(wlp::Shape::CIRCLE);
        } else if (classId == -3) {
            mref.size(Size(0.5,.5,0.)).color(Color::Black)
            .shape(wlp::Shape::CIRCLE);
        }
        else {
            // mref.color(getColorFromPosition(position.x(), position.y(), targetResolution));
            mref.color(getColorFromClassId(classId));
        }

        m = mref.id();
            }
    void del(){
        wlp::Mark &mref = dbVis->getMark(m);
        mref.del();
    }
    void updatePosition(Vec2d &position) {
        wlp::Mark &mref = dbVis->getMark(m);
        mref.position(position.x(), position.y(), 0);
    }
    uint m;
    wlp::Fatum *dbvis;
};

/**
 * @brief The PolyM struct helps drawing borders given by the Voronoi tesselation
 */
struct PolyM {
    PolyM(wlp::Fatum * dbvis, const vector<Vec2d> &poly):dbvis(dbvis) {
        dbvis->defaultMark().color(Color::Black).size(Size(0.01,0.01,0.01)).shape(wlp::Shape::CIRCLE);
        for (unsigned i=0; i<poly.size(); ++i)
            vertex.push_back(dbvis->addMark().position(poly[i].x(), poly[i].y(), 0));
        for (unsigned i=0; i<vertex.size(); ++i) {
            Connection con = dbvis->addConnection(vertex[i], vertex[(i+1)%poly.size()]);
            con.sourceColor(Color::Black);
            con.targetColor(Color::Black);
            con.arrow(true);
            lines.push_back(con);
        }
    }
    void del() {
        for (unsigned i=0; i<vertex.size(); ++i)
            vertex[i].del();
    }
    vector<wlp::Mark> vertex;
    vector<wlp::Connection> lines;
    wlp::Fatum * dbvis;
};

/**
 * @brief outputPositions : write positions in a csv type file
 * @param pos : positions of the vertices
 * @param vtt : vector to tell which vertex to write (true) or not (false)
 * @param header : information to write out before writing positions (for filter split)
 * @param filename : filename/path where to write the positions
 */
void outputPositions(const std::vector<Vec2d> &pos, const std::vector<int> &vtt, std::string filename) {
    std::ofstream f;
    f.open(filename.c_str(), ios::out);
    for (unsigned i = 0; i < pos.size(); ++i) {
        if (vtt[i] > -1) {
            // f << pos[i].x() << "," << pos[i].y() << std::endl;
            f << i << "," << pos[i].x() << "," << pos[i].y() << std::endl;
        }
    }
    f.close();
}

/**
 * @brief normalizePositions : normalize given positions into a resolutionW x resolutionH space
 * @param pos : positions to normalize
 * @param resolutionW : width of the space to normalize to
 * @param resolutionH : height of the space to normalize to
 * @param keepBothRatio : preserve the shape or not of the positions
 */
void normalizePositions(std::vector<Vec2d> &pos, int resolutionW, int resolutionH, bool keepBothRatio) {
    wlp::Vec2d minPos(pos[0]);
    wlp::Vec2d maxPos(pos[0]);
    for (unsigned i = 0; i < pos.size(); ++i) {
        Vec2d n = pos[i];
        minPos = wlp::minVector(minPos, n);
        maxPos = wlp::maxVector(maxPos, n);
    }
    Vec2d offset = minPos;
    Vec2d delta = maxPos - offset;
    if (!keepBothRatio) {
        double maxDim = delta.x() > delta.y() ? delta.x() : delta.y();
        delta = Vec2d(maxDim, maxDim);
    }
    Vec2d factor(resolutionW-2, resolutionH-2);
    for (unsigned i = 0; i < pos.size(); ++i) {
        pos[i] = (pos[i] - offset);
        pos[i].setX(pos[i].x() / delta.x() * factor.x() + 1);
        pos[i].setY(pos[i].y() / delta.y() * factor.y() + 1);
    }
}

/**
 * @brief centerPositions : center vertex positions around (0, 0)
 * @param pos : vertex positions to center
 * @param resolutionW : width of the centered space
 * @param resolutionH : height of the centered space
 */
void centerPositions(std::vector<Vec2d> &pos, int resolutionW, int resolutionH) {
    double midW = resolutionW * 1.0 / 2;
    double midH = resolutionH * 1.0 / 2;
    Vec2d offset(midW, midH);
    for (unsigned i = 0; i < pos.size(); ++i) {
        pos[i] -= offset;
    }
}

/**
 * @brief fixOverlap : avoid overlapping vertices with random (should be changed in a deterministic way)
 * @param pos : positions to check for overlapping
 * @param radius : radius of the random wobble applied to the overlapping vertices
 */
void fixOverlap(std::vector<Vec2d> &pos, std::vector<int> classIds, double radius = 1e-3) {
    IncrementalDelaunay incr;
    for (unsigned i = 0; i < pos.size(); ++i) {
        try {
            incr.addPoint(pos[i]);
        } catch (...){
            double rPosX = ((double)rand() / RAND_MAX)*2*radius - radius;
            double rPosY = ((double)rand() / RAND_MAX)*2*radius - radius;
            pos[i] += Vec2d(rPosX, rPosY);
            --i;
        }
    }
}

/**
 * @brief borderPoints : build a border of dummy vertices around a resolutionW x resolutionH rectangle
 * @param resolutionW : width of the rectangle
 * @param resolutionH : height of the rectangle
 * @return : vector of dummy vertices
 */
vector<Vec2d> borderPoints(int resolutionW, int resolutionH) {
    vector<Vec2d> result( 4 + resolutionW * 2 + resolutionH * 2);
    double topCoordinate = resolutionH * 1.0 / 2 + 0.5;
    double bottomCoordinate = -topCoordinate;
    double rightCoordinate = resolutionW * 1.0 / 2 + 0.5;
    double leftCoordinate = -rightCoordinate;
    int index = 0;
    result[index++] = Vec2d(leftCoordinate, topCoordinate);
    result[index++] = Vec2d(leftCoordinate, bottomCoordinate);
    result[index++] = Vec2d(rightCoordinate, topCoordinate);
    result[index++] = Vec2d(rightCoordinate, bottomCoordinate);

    for (int x = 1; x <= resolutionW; ++x) {
        result[index++] = Vec2d(leftCoordinate + x, topCoordinate);
        result[index++] = Vec2d(leftCoordinate + x, bottomCoordinate);
    }
    for (int y = 1; y <= resolutionH; ++y) {
        result[index++] = Vec2d(leftCoordinate , bottomCoordinate + y);
        result[index++] = Vec2d(rightCoordinate, bottomCoordinate + y);
    }

    return result;
}

std::vector<int> fixedNodes;
int fixedNodesCount = 0;
/**
 * @brief buildBorder : calls borderPoints and keep track of those dummy vertex IDs into fixedNodes
 * @param nodePos : positions of samples to mix dummy vertices into
 * @param resolutionW : width of the target resolution
 * @param resolutionH : height of the target resolution
 */
void buildBorder(std::vector<Vec2d> &nodePos, std::vector<int> &indexesToTake, int resolutionW, int resolutionH) {
    std::vector<Vec2d> bp = borderPoints(resolutionW, resolutionH);
    int len = bp.size() + indexesToTake.size();
    fixedNodes = std::vector<int>(len, 0);
    // for (int i = 0; i < len; ++i) {
    //     fixedNodes[i] = 0;
    // }
    int base = indexesToTake.size();
    for (int ip = 0; ip < bp.size(); ++ip ) {
        Vec2d &p = bp[ip];
        fixedNodes[base + ip] = 1;
        fixedNodesCount += 1;
        nodePos.push_back(p);
        indexesToTake.push_back(-1);
    }
}

/**
 * @brief normalizeActiveSampleNodes : normalize not fixed vertices into the remaning square space. Checked with fixedNodes.
 * @param nodePos : vertices, fixed and/or not, to normalize
 * @param currentResolution : width/height of the remaining empty space.
 */
void normalizeActiveSampleNodes(std::vector<Vec2d> &nodePos, int currentResolution) {
    double halfResolution = currentResolution * 1.0 / 2 + 0.5;
    for (unsigned i = 0; i < nodePos.size(); ++i) {
        if (fixedNodes[i] == 0) {
            Vec2d newPos = nodePos[i] * (halfResolution - 1) / halfResolution;
            nodePos[i] = newPos;
        }
    }
}

std::vector<PolyM> polygonShapes;
/**
 * @brief drawPolygon : draws polygon given a vector of vertices and two colors. Directly added into polygonShapes vector.
 */
void drawPolygon(std::vector<Vec2d> &polygon, wlp::Color vertexColor, wlp::Color edgeColor) {
    polygonShapes.push_back(PolyM(dbVis, polygon));
}

/**
 * @brief clearPolygons : remove every drawn polygon by drawPolygon function
 */
void clearPolygons() {
    for (auto &pm : polygonShapes) {
        pm.del();
    }
    polygonShapes.clear();
}

std::vector<SampleMark> sampleMarks;
/**
 * @brief initializeSampleMarks : prepare an array of SampleMark to draw each sample, colored by their GT and checking if they are well predicted or not.
 * @param points : position of every vertice to draw
 * @param classIds : sample GTs
 * @param classPreds : sample predictions
 */
void initializeSampleMarks(std::vector<Vec2d> &points, std::vector<int> &indexesToTake, int targetResolution) {
    for (auto &m : sampleMarks) {
        m.del();
    }

    sampleMarks.clear();
    for (unsigned i = 0; i < points.size(); ++i) {
        sampleMarks.push_back(SampleMark(dbVis, points[i], indexesToTake[i], targetResolution)); //border and dummy nodes
    }
}

/**
 * @brief updateSampleMarks : update SampleMark positions
 * @param points : new positions of all the vertices
 */
void updateSampleMarks(std::vector<Vec2d> &points) {
    for (unsigned i = 0; i < sampleMarks.size(); ++i) {
        sampleMarks[i].updatePosition(points[i]);
    }
}

int currentResolution = -1;
/**
 * @brief initializeRoutine : prepare algorithm by placing vertices in 2d space along with border
 * @param pos : initial vertex positions
 * @param resolutionW : width of target resolution
 * @param resolutionH : height of target resolution
 */
void initializeRoutine(std::vector<Vec2d> &pos, std::vector<int> &indexesToTake, int resolutionW, int resolutionH) {
    fixOverlap(pos, indexesToTake);
    buildBorder(pos, indexesToTake, resolutionW, resolutionH);
    initializeSampleMarks(pos, indexesToTake, resolutionW);
    currentResolution = resolutionW;
}


/**
 * @brief completeWithDummyPoints : if there is no enough samples to fill a square, create dummy nodes and return them
 * @param numberOfActiveSamples : number of "true" samples
 * @param targetResolution : width/height of the targeted square space;
 * @return : a vector with the dummy sample positions
 */
std::vector<Vec2d> completeWithDummyPoints(int numberOfActiveSamples, int targetResolution) {
    int delta = targetResolution * targetResolution - numberOfActiveSamples;
    std::vector<Vec2d> result(delta);
    for (int dummyIndex = 0; dummyIndex < delta; ++dummyIndex) {
        int toX = dummyIndex % 4 < 2 ? -1 : 1;
        int toY = dummyIndex % 4 == 0  || dummyIndex % 4 == 3 ? 1 : -1;
        int dX = 0, dY = 0;
        if (toX < 0 && toY < 0) {
            dX = 1;
        }
        if (toX > 0 && toY > 0) {
            dX = -1;
        }
        if (toX > 0 && toY < 0) {
            dY = 1;
        }
        if (toX < 0 && toY > 0) {
            dY = -1;
        }
        // double rPosX = ((double)rand() / RAND_MAX)*targetResolution;
        // double rPosY = ((double)rand() / RAND_MAX)*targetResolution;
        result[dummyIndex] = Vec2d(toX * targetResolution / 2, toY * targetResolution / 2) + Vec2d(dX * floor(dummyIndex / 4), dY * floor(dummyIndex/4));
    }

    return result;
}



typedef std::vector<std::vector<bool>> virtualGrid;
std::vector<Vec2d> fillGapWithVirtualGrid(std::vector<Vec2d> &nodePos, int targetResolution, int internalResolution) {
    virtualGrid vg(internalResolution);
    std::vector<Vec2d> result;
    for (int i = 0; i < internalResolution; ++i) {
        vg[i] = std::vector<bool>(internalResolution, false);
    }

    for (int i = 0; i < nodePos.size(); ++i) {
        int x = std::floor((nodePos[i].x() / (double)targetResolution + 0.5) * (double)internalResolution);
        int y = std::floor((nodePos[i].y() / (double)targetResolution + 0.5) * (double)internalResolution);
        vg[x][y] = true;
    }

    for (int i = 0; i < internalResolution; ++i) {
        for (int j = 0; j < internalResolution; ++j) {
            if (!vg[i][j]) {
                double x = ((((double)i) /(double)internalResolution) - 0.5) *(double)targetResolution;
                double y = ((((double)j) /(double)internalResolution) - 0.5) *(double)targetResolution;
                result.push_back(Vec2d(x, y)); 
            }
        }
    }

    return result;
}

double triangleArea(Vec2d a, Vec2d b, Vec2d c) {
    double d1 = (a - b).length();
    double d2 = (a - c).length();
    double d3 = (b - c).length();
    double p = (d1 + d2 + d3)/2.;
    return sqrt(p * (p - d1) * (p - d2) * (p - d3));
}


const float fillFactor = .025;
// higher factor = more points
std::vector<Vec2d> fillGapWithDelaunay(std::vector<Vec2d> &nodePos, int targetResolution, int internalResolution) {
    std::vector<Vec2d> result;
    IncrementalDelaunay incr;
    auto ids = incr.getIds();
    
    for (int i = 0; i < nodePos.size(); ++i) {
        Vec2d &p = nodePos[i];
        incr.addPoint(p, i);
    }

    std::vector<IncrementalDelaunay::Triangle> tris = incr.getTriangles();
    for (int i = 0; i < tris.size(); ++i) {
        IncrementalDelaunay::Triangle t = tris[i];
        Vec2d a = nodePos[ids[t.n1]];
        Vec2d b = nodePos[ids[t.n2]];
        Vec2d c = nodePos[ids[t.n3]];
        double area = triangleArea(a, b, c);
        cerr << area << endl;
        Vec2d center = (a + b + c)/3.;
        for (int j = 1; j <= area/targetResolution/fillFactor; ++j) {
            double rPosX = ((double)rand() / RAND_MAX)*2 - 1;
            double rPosY = ((double)rand() / RAND_MAX)*2 - 1;
            Vec2d v = Vec2d(center.x() + rPosX, center.y() + rPosY);
            result.push_back(v);
        }
    }

    return result;
}

const float distanceCut = 5.;
std::vector<Vec2d> fillGapWithDelaunayOnDistances(std::vector<Vec2d> &nodePos, int targetResolution, int internalResolution) {
    std::vector<Vec2d> result;
    IncrementalDelaunay incr;
    auto ids = incr.getIds();
    
    for (int i = 0; i < nodePos.size(); ++i) {
        Vec2d &p = nodePos[i];
        incr.addPoint(p, i);
    }

    std::vector<IncrementalDelaunay::Triangle> tris = incr.getTriangles();
    for (int i = 0; i < tris.size(); ++i) {
        IncrementalDelaunay::Triangle t = tris[i];
        Vec2d a = nodePos[ids[t.n1]];
        Vec2d b = nodePos[ids[t.n2]];
        Vec2d c = nodePos[ids[t.n3]];
        double area = triangleArea(a, b, c);
        float ab = (a-b).length();
        float ac = (a-c).length();
        float bc = (b-c).length();
        double rPosX = ((double)rand() / RAND_MAX)*.25 - .125;
        double rPosY = ((double)rand() / RAND_MAX)*.25 - .125;
        Vec2d rp(rPosX, rPosY);
        if (ab > distanceCut)
            result.push_back(a+(a-b)/2. + rp);
        if (ac > distanceCut)
            result.push_back(a+(a-c)/2. + rp);
        if (bc > distanceCut)
            result.push_back(b+(b-c)/2. + rp);
    }
    return result;
}

std::vector<Vec2d> fillGap(std::vector<Vec2d> &nodePos, int targetResolution, int internalResolution) {
    internalResolution = targetResolution * 0.9;
    // return fillGapWithDelaunayOnDistances(nodePos, targetResolution, internalResolution);
    // return fillGapWithDelaunay(nodePos, targetResolution, internalResolution);
    // return fillGapWithVirtualGrid(nodePos, targetResolution, internalResolution);
    return std::vector<Vec2d>();
}

// /**
//  * @brief prepareIteration : do the computation to filter and fill with dummy samples
//  * @param nodePos : original vertex positions
//  */
// void prepareIteration(std::vector<Vec2d> &nodePos,
//                       std::vector<int> &indexesToTake) {
//     nodePos.clear();
//     indexesToTake.clear();
//     int actualCount = 0;
//     actualCount = nodePosOrig.size();
//     indexesToTake = std::vector<int>(actualCount);
//     std::iota(indexesToTake.begin(), indexesToTake.end(), 0);
//     // for (int i = 0; i < actualCount; ++i){
//     //     indexesToTake[i] = i;
//     // }

//     nodePos = std::vector<Vec2d>(actualCount);
//     for (int i = 0; i < actualCount; ++i) {
//         nodePos[i] = nodePosOrig[indexesToTake[i]];
//     }

//     int targetResolution = ceil(sqrt(actualCount));
//     normalizePositions(nodePos, targetResolution, targetResolution, true);
//     centerPositions(nodePos, targetResolution, targetResolution);
//     std::vector<Vec2d> dummies = completeWithDummyPoints(actualCount, targetResolution);
//     std::vector<Vec2d> spaceFillers = fillGap(nodePos, targetResolution, targetResolution); // magic number here
//     int addedFillers = spaceFillers.length();
//     for (unsigned i = 0; i < dummies.size(); ++i) {
//         nodePos.push_back(dummies[i]);
//         indexesToTake.push_back(-2);
//     }
//     for (unsigned i = 0; i < spaceFillers.size(); ++i) {
//         nodePos.push_back(spaceFillers[i]);
//         indexesToTake.push_back(-3);
//     }

//     targetResolution += addedFillers
//     initializeRoutine(nodePos, indexesToTake, targetResolution, targetResolution);
// }

vector<int> readGT(string gtFile) {
    std::ifstream in(gtFile.c_str());
    char line[255];
    int index = 0;
    vector<int> result(0);
    while (!in.eof()) {
        in.getline(line, 255);
        string lines(line);
        // vector<string> token;
        // tokenize(lines, token, ",");
        // if (token.size() >= 0) {
            int gt = atoi(line);
            // weak
            result.push_back(gt);
            ++index;
        // }
    }

    return result;
}

/**
 * @brief prepareIteration : do the computation to filter and fill with dummy samples
 * @param nodePos : original vertex positions
 */
void prepareIteration(std::vector<Vec2d> &nodePos,
                      std::vector<int> &indexesToTake) {
    nodePos.clear();
    indexesToTake.clear();
    int actualCount = 0;
    actualCount = nodePosOrig.size();
    indexesToTake = std::vector<int>(actualCount);
    // std::iota(indexesToTake.begin(), indexesToTake.end(), 0);
    vector<int> gts = readGT("labels-swiss.csv");
    for (int i = 0; i < actualCount; ++i){
        // indexesToTake[i] = i;
        indexesToTake[i] = gts[i];
        // cerr << gts[i] << endl;
    }


    nodePos = std::vector<Vec2d>(actualCount);
    for (int i = 0; i < actualCount; ++i) {
        // nodePos[i] = nodePosOrig[indexesToTake[i]];
        nodePos[i] = nodePosOrig[i];
    }

    int targetResolution = ceil(sqrt(actualCount));
    normalizePositions(nodePos, targetResolution, targetResolution, true);
    centerPositions(nodePos, targetResolution, targetResolution);
    std::vector<Vec2d> spaceFillers = fillGap(nodePos, targetResolution, targetResolution); // magic number here
    int addedFillers = spaceFillers.size();
    cerr << addedFillers << endl;
    actualCount += addedFillers;
    targetResolution = ceil(sqrt(actualCount));
    std::vector<Vec2d> dummies = completeWithDummyPoints(actualCount, targetResolution);
    for (unsigned i = 0; i < dummies.size(); ++i) {
        nodePos.push_back(dummies[i]);
        indexesToTake.push_back(-2);
    }
    for (unsigned i = 0; i < spaceFillers.size(); ++i) {
        nodePos.push_back(spaceFillers[i]);
        indexesToTake.push_back(-3);
    }

    normalizePositions(nodePos, targetResolution, targetResolution, true);
    centerPositions(nodePos, targetResolution, targetResolution);
    initializeRoutine(nodePos, indexesToTake, targetResolution, targetResolution);
}



/**
 * @brief arrangeIter : core algorithm
 *              block a : - compute voronoy tesselation and store each face in 2d vectors
 *                        - for each face :
 *                              - draw it
 *                              - check if its centroid is inside the targeted space
 *                                  - if not : do not move the concerned vertex for this iteration
 *                                  - if yes and if concerned vertex is not fixed : move the concerned vertex to the centroid and update delta values
 *
 *              block b : - check if there was two centroids outside of targeted space
 *                          - if yes : scale down active vertices into targeted space for retry
 *
 *              block c : - if all vertices could be moved to their concerned centroids and movements were somewhat minimal:
 *                          - place every active vertex into IncrementalDelaunay structure
 *                          - compute places to snap the most outward vertices to an inside-border in the free space
 *                          - for each of those places :
 *                             - move the nearest active vertex to that place, make it fixed
 *                          - reduce free space resolution and go back to block a until there is no free space left
 * @param pos : each vertex positions (true, dummy, fixed)
 * @param iter : current iteration number (logging purpose)
 * @return : average movement (delta) among every active vertices
 */
int subIter = 0;
Vec2d arrangeIter(std::vector<Vec2d> &pos, long iter) {
    Vec2d deltaSum = Vec2d(0, 0);
    Vec2d deltaMax = Vec2d(0, 0);
    Voronoy vor;
    std::vector<Vec2d> vPoints;                     // contains a list of points used to draw faces
    std::vector<std::vector<unsigned int>> vFaces;  // ordered vPoints indexes for each face computed
    std::vector<unsigned int> vClippingPoints;      // list of vPoints indexes that are out of bounds
    vor.voronoyFromPoints(pos, vPoints, vFaces, vClippingPoints);

    double xymin = -(currentResolution-0.5) * 1.0 / 2;
    double xymax = -xymin;
    wlp::Rectd border = wlp::Rectd(xymin, xymin, xymax, xymax);

    clearPolygons();
    int gridcontact = 0;

    // BLOCK A
    for (unsigned fId = 0; fId < vFaces.size(); ++fId) {
        if (fixedNodes[fId] == 0) {
            std::vector<Vec2d> polygon(vFaces[fId].size());
            for (unsigned pId = 0; pId < vFaces[fId].size(); ++pId) {
                polygon[pId] = vPoints[vFaces[fId][pId]];
            }
            drawPolygon(polygon, wlp::Color::Orange, wlp::Color::Green);
            wlp::Poly2d p = wlp::Poly2d(polygon);
            Vec2d center = p.getBarycenter();

            if (!border.isInside(center)) {
                //cout << "[Warning] move out of grid border" << endl;
                center = pos[fId];
                ++gridcontact;
            }

            if (fixedNodes[fId] == 0) {
                Vec2d delta = (center - pos[fId]);
                Vec2d clamp = delta / delta.length() * maxRange;
                if (delta.length() < maxRange) {
                    pos[fId] = center;
                    delta = Vec2d(abs(delta.x()), abs(delta.y()));
                } else {
                    pos[fId] = pos[fId] + clamp;
                    delta = Vec2d(abs(clamp.x()), abs(clamp.y()));
                }
                deltaSum += delta;
                deltaMax = maxVector(deltaMax, delta);
            }
        }
    }

    cerr << deltaSum.length() << endl;
    
    // BLOCK B
    if (gridcontact > 0) {
        for (unsigned nId = 0; nId < pos.size(); ++nId) {
            if (fixedNodes[nId] == 0) {
                pos[nId] *= currentResolution /(currentResolution+1.);
            }
        }
    }


    // deltaSum /=  vFaces.size()-fixedNodes.size();
    // std::cerr << iter << " (" << vFaces.size() - fixedNodes.size() <<  " | " << fixedNodes.size() << ") - dSum : " << deltaSum << " dMax:" << deltaMax<< std::endl;

    //BLOCK C
    // if (gridcontact == 0 && deltaMax.length() < 2e-2 && deltaSum.length() < 1e-2 && currentResolution > 1) {
    if (useIterativeBorder) {
        bool condition = false;
        if (useFixedIterationNumber) {
            ++subIter;
            int numActive = vFaces.size() - fixedNodesCount;
            condition = subIter >= sqrt(numActive)/iterationInverseFactor && currentResolution >= 1;
        } else {
            condition = gridcontact == 0 && deltaMax.length() < sensibility && currentResolution >= 1;
        }



        if (condition) { 
            if (useFixedIterationNumber) {
                subIter = 0;
            } else if (iter%snapCheckFreq == 0){

            //use delauny triangulation to find nearest neighbour in O(1)
            IncrementalDelaunay delaunay;
            auto ids = delaunay.getIds();

            for (unsigned nId = 0; nId < pos.size(); ++nId) {
                if (fixedNodes[nId] == 0) {
                    delaunay.addPoint(pos[nId], nId);
                }
            }
            std::vector<Vec2d> snapPoints = borderPoints(currentResolution - 2, currentResolution - 2);
            std::vector<int> candidates(snapPoints.size(), -1);
            std::vector<int> tempFixedNodes = fixedNodes;
            // for (const auto &p: snapPoints) {
            bool success = true;
            int candidateCount = 0;
            for (int i = 0; i < snapPoints.size(); ++i) {
                auto p = snapPoints[i];
                node n = delaunay.getClosestNode(p);
                while (fixedNodes[ids[n]] == 1 || tempFixedNodes[ids[n]] == 1) {
                    delaunay.getClosestPoint(p);
                    delaunay.delNode(n);
                    n = delaunay.getClosestNode(p);
                }
                candidates[i] = ids[n];
                tempFixedNodes[ids[n]] = 1;
                double dist = (pos[ids[n]] - p).length();
                double currentRatio = i/((double)candidates.size());
                if (dist > maxRangeSnap && currentRatio < forceSnapRatio) {
                    success = false;
                    candidateCount = i;
                    i = snapPoints.size();
                }
            }
            double ratio = ((double) candidateCount)/((double)candidates.size());
            cerr << candidateCount << "/" << candidates.size() << " : " << ratio << endl;
            if (success) {
                for (int i = 0; i < candidates.size(); ++i) {
                    int ind = candidates[i];
                    pos[ind] = snapPoints[i];
                    fixedNodes[ind] = 1;
                    fixedNodesCount += 1;
                }
                
                currentResolution -= 2;
            }

            if (currentResolution == 1) { // only the 0, 0 node to snap
                node n = delaunay.getClosestNode(Vec2d(0, 0));
                pos[ids[n]] = Vec2d(0, 0);
                currentResolution = 0;
            }
            }
        }
    }



    return deltaMax;
}

void removeSpaceFillers(std::vector<Vec2d> &nodes, std::vector<int> &nodeIds) {
    std::vector<Vec2d> newNodes();
    std::vector<int> newIds();
    int size = nodeIds.size();
    for (int i = 0; i < size; ++i) {
        if (nodeIds[i] == -3) {
            nodes.erase(nodes.begin() + i);
            nodeIds.erase(nodeIds.begin() + i);
            fixedNodes.erase(fixedNodes.begin() + i);
            i -= 1;
            size -= 1;
        }
    }
    
}




Vec2d lastDelta = Vec2d(1,1);
std::vector<Vec2d> nodePos;
std::vector<int> nodesToTake;
long iter = 0;
bool fini = false;
bool hasStabilizedOnce = false;
bool started = false;
static void glut_idle_callback(void) {
    if(dbVis->needsRendering()) {
        glutPostRedisplay();
        return;
    }
    if (!started) {
        return;
    }
    // if (iter == 0) {
    //     std::this_thread::sleep_for(std::chrono::milliseconds(4000));
    // }
    // if (hasStabilizedOnce) {
    //     std::this_thread::sleep_for(std::chrono::milliseconds(4000));
    //     hasStabilizedOnce = false;
    // }
    if (!fini) {
        if (currentResolution >= 0) {
            lastDelta = arrangeIter(nodePos, ++iter);
            if (lastDelta.x() < gapSensibility && lastDelta.y() < gapSensibility && !hasStabilizedOnce) {
                removeSpaceFillers(nodePos, nodesToTake);
                hasStabilizedOnce = true;
                useFixedIterationNumber = false;
                useIterativeBorder = true;
            }
            // outputPositions(nodePos, nodesToTake, std::to_string(iter));
            if (iter % 100 == 0) {
            }
            fini = currentResolution < 1;
            if (fini) {
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed_seconds = end - beginC;
                std::cerr << "End ! Took " << elapsed_seconds.count() << "s." << std::endl;

            }
        }
        if (currentResolution < 0){

        } else {
            updateSampleMarks(nodePos);
        }
    } else {
    }
    dbVis->swap();
    // garder les anciennes marks
    //dbVis->animate(wlp::AnimParameters(100));
    glutPostRedisplay();
    return;
}


static void glut_reshape_callback(int width, int height) {
    dbVis->getCamera().setViewport(Vec4i(0,0,width,height));
    WIDTH = width;
    HEIGHT = height;
    glutPostRedisplay();
}
static void glut_draw_callback(void) {
    dbVis->renderFrame();
    glutSwapBuffers();
    if(dbVis->needsRendering()) glutPostRedisplay();
}
static void glut_special(int special, int , int ) {
    switch (special) {
    default:
        break;
    }
    glutPostRedisplay();
}
static void keyboard(const unsigned char key, const int , const int ) {
    switch(int(key)) {
    case 'q':
        outputPositions(nodePos, nodesToTake, outputFilename);
        std::cout << "EXIT ... " << std::endl;
        exit(EXIT_SUCCESS);
    case 's':
        started = true;
    case 'n':
        fini = false;
    default:
        return;
    }
    glutPostRedisplay();
}
static void motionFunc(int x, int y) {
    //if (dbVis->needsRendering()) return; //do not treat events until the graph is re rendered (else too many events)
    if (x < 0 || x > WIDTH || y < 0 || y>HEIGHT) return;
    if (lastX < 0 || lastX>WIDTH || lastY < 0 || lastY>HEIGHT) return;
    if (mouseState == GLUT_LEFT_BUTTON) {
        Vec2f current = Vec2f(x, y);
        Vec3f currentModel = dbVis->getCamera().windowToModel(current);
        Vec3f previousModel = dbVis->getCamera().windowToModel(Vec2f(lastX, lastY));
        Vec3f move = currentModel-previousModel;
        dbVis->getCamera().moveModel(move[0],move[1]);
        dbVis->getCamera().swap();
        glutPostRedisplay();
        lastX = x;
        lastY = y;
    }
}
static void mouseFunc(int button, int state, int x, int y) {
    if (x < 0 || x > WIDTH || y < 0 || y>HEIGHT) return;
    mouseState = button;
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            lastX = x;
            lastY = y;
        }
        // wheel up
    } else if (button == 3) {
        dbVis->getCamera().zoom(1.05, Vec2f(x,y));
        dbVis->getCamera().swap();
        if(dbVis->needsRendering()) glutPostRedisplay();
    } else if (button == 4) {
        dbVis->getCamera().zoom(0.95, Vec2f(x,y));
        dbVis->getCamera().swap();
        if(dbVis->needsRendering()) glutPostRedisplay();
    } else if (button == GLUT_RIGHT_BUTTON) {
        lastX = x;
        lastY = y;
        if(state == GLUT_DOWN) {
            motionFunc(x,y);
        }
    }
}
/**
 * @brief initVis : initialize Fatum parameters
 */
void initVis() {
    dbVis = new wlp::Fatum();
    dbVis->initGL();
    dbVis->getCamera().setViewport(Vec4i(0,0,100,100));
    dbVis->layerOn(MARKS | TEXT | CONNECTIONS);
}

/**
 * @brief readNodePositions : read csv type file to place every vertex in a 2d space
 *      file format : 1 vertex per line : <x>,<y>
 * @param filename : csv filename/path
 * @param positions : pre-initialized vector to contains positions of read vertices
 */
void readNodePositions(const std::string filename, std::vector<Vec2d> &positions, int sampleCount) {
    io::CSVReader<2, io::trim_chars<'\t'>, io::no_quote_escape<','>>
        inN(filename);
    int index = 0;
    double posX, posY;
    while(inN.read_row(posX, posY)) {
        if (index < sampleCount) {
            positions[index] = Vec2d(posX, posY);
        }
        ++index;
    }
}

int main(int argc, char* argv[]) {
    std::cerr << "Usage : vrgrid sample_positions sample_count placement_result" << std::endl;
    glutInit(&argc, argv);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH | GLUT_MULTISAMPLE);
    if ((glutCreateWindow("Wulip Glut Viewer")) == GL_FALSE) {
        std::cerr << "Unable to create a OpenGl Glut window" << std::endl;
        exit(EXIT_FAILURE);
    }
    glutIdleFunc   (glut_idle_callback      );
    glutReshapeFunc(glut_reshape_callback   );
    glutDisplayFunc(glut_draw_callback      );
    glutSpecialFunc(glut_special);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);

    initVis();

    std::string positionsFilename = argv[1];
    int sampleCount = atoi(argv[2]);
    outputFilename = argv[3];
    nodePosOrig = std::vector<Vec2d>(sampleCount);
    
    readNodePositions(positionsFilename, nodePosOrig, sampleCount);
    beginC = std::chrono::steady_clock::now();
    prepareIteration(nodePos, nodesToTake);

    dbVis->swap();
    dbVis->center();
    dbVis->getCamera().swap();
    glutPostRedisplay();
    glutMainLoop();
    return 0;
}
