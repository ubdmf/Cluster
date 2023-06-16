#include <iostream>
#include <vector>
#include <cmath>
#include<deque>
#include<fstream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <sstream>
#include <algorithm>
#include<stack>
#include<map>

using namespace std;
std::ofstream labels("/Users/xinyi.wang2/CLionProjects/untitled/label2.csv");
std::ofstream polygen_test2("/Users/xinyi.wang2/CLionProjects/untitled/polygen2.csv");
std::ofstream final_polygen2("/Users/xinyi.wang2/CLionProjects/untitled/final_polygen2.csv");
std::ofstream orig_point2("/Users/xinyi.wang2/CLionProjects/untitled/ori_point2.csv");
std::ofstream radius_test2("/Users/xinyi.wang2/CLionProjects/untitled/radius2.csv");
std::ofstream k1_test("/Users/xinyi.wang2/CLionProjects/untitled/k1beforesort.csv");
std::ofstream k2_test("/Users/xinyi.wang2/CLionProjects/untitled/k2beforesort.csv");
std::ofstream kd1_test("/Users/xinyi.wang2/CLionProjects/untitled/kd1test.csv");
class Point {
public:
    float x = 0;
    float y = 0;
    int cluster = 0;
    int type = 0; //  3 = core, 2 = boundry 1 = noise
};

class DBSCAN {
public:
    std::vector<Point> all_pt_;
    std::vector<Point> core_pt;
    std::vector<int> label_;
    std::vector<vector<float>> distance_mat_;
    int min_pts_;
    float radius_;

    DBSCAN() = default;

    DBSCAN(std::vector<Point> all_pts, int min_pts, float radius) : all_pt_(all_pts), min_pts_(min_pts),
                                                                    radius_(radius) {
    };

    DBSCAN(std::vector<Point> all_pts, int min_pts) : all_pt_(all_pts), min_pts_(min_pts) {
        radius_ = CalcEpsilon2(all_pts, min_pts);
        std::cout << "radius " << radius_ << std::endl;
    }

    void JoinEdgePoint();

    ~DBSCAN() {};

    void GetCluster();

    bool ExpandCluster(int pt_idx, int &cluster_idx);

    float CalcDistance(Point &pt1, Point &pt2);

    float eval();

    void viz();

    //可以使用肘点检测方法来得出一个合适的epsilon值。
    //计算每个点和它的k个最近的邻居之间的平均距离，其中k=我们选择的MinPts。然后，我们将这些平均k距离按升序绘制在k距离图上。
    //epsilon的最佳值是具有最大曲率或弯曲的点，即在最大的斜率处。
    float CalcEpsilon2(std::vector<Point> &all_pt, int MinPts);


};


void DBSCAN::JoinEdgePoint() {
    for (int i = 0; i < all_pt_.size(); i++) {
        if (label_[i] == -1) {
            for (int j = 0; j < all_pt_.size(); j++) {
                if (label_[j] > 0) {
                    if (distance_mat_[i][j] <= radius_) {
                        label_[i] = label_[j];
                    }
                }
            }
        }
    }

}

float DBSCAN::eval() {
    int count = 0;
    for (int i = 0; i < all_pt_.size(); i++) {
        if (label_[i] == -1) {
            count++;
        }
    }
    return count;
}

//计算每个点附近相邻10个点的平均距离
//n2耗时
//map 里存取的是每个点距离其他所有点的range，可以选取minpts个
float DBSCAN::CalcEpsilon2(std::vector<Point> &all_pt, int MinPts) {
    MinPts = 5;
    std::vector<float> distance_mat_ave;
    vector<vector<std::pair<int, float>>> distance_mat_2;
    for (int i = 0; i < all_pt.size(); i++) {
        distance_mat_ave.resize(all_pt.size());
        distance_mat_2.resize(all_pt.size());
        for (int j = 0; j < all_pt.size(); j++) {
            if (i == j) { continue; }
            float x1 = all_pt[i].x;
            float y1 = all_pt[i].y;
            float x = all_pt[j].x;
            float y = all_pt[j].y;
            float diff_range = std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
            distance_mat_2[i].push_back(std::make_pair(j, diff_range));
        }
        sort(distance_mat_2[i].begin(), distance_mat_2[i].end(),
             [](const pair<int, float> &a, pair<int, float> &b) { return a.second < b.second; });
        distance_mat_2[i].resize(MinPts);
        float sum = 0;
        for (int n = 0; n < MinPts; n++) {
            sum += distance_mat_2[i][n].second;
            std::cout << " DEBUG " << distance_mat_2[i][n].second << std::endl;
            std::cout << "DEBUG " << i << " DEBUG  Sum" << sum << std::endl;
        }
        float ave = sum / MinPts;
        distance_mat_ave[i] = ave;
    }

    sort(distance_mat_ave.begin(), distance_mat_ave.end(), [](float a, float b) { return a > b; });
    //两两之间求斜率，找到delta k 变化最大的点
    int kne = 0;
    int max = -10000;
    std::vector<float> div_k;

//    std::vector<float> div_k;
    for (int i = 5; i < distance_mat_ave.size()-5; i++) {
        float kd2 = std::atan((distance_mat_ave[i + 2] - distance_mat_ave[i + 1])) -
                             std::atan((distance_mat_ave[i + 1] - distance_mat_ave[i]))/
                    std::sqrt((distance_mat_ave[i + 1] - distance_mat_ave[i])*(distance_mat_ave[i + 1] - distance_mat_ave[i])+1);
        float kd3 = (distance_mat_ave[i + 5] + distance_mat_ave[i -5] - 2* distance_mat_ave[i])/(2*5);
        k2_test << kd3 << std::endl;

        k1_test << kd2 << std::endl;
        radius_test2 << distance_mat_ave[i] << std::endl;
        float k_d1 = (distance_mat_ave[i - 5] - distance_mat_ave[i+5]) / 2;
        kd1_test << k_d1 << std::endl;
        div_k.push_back(k_d1);
    }
//    std::map<float,int> div2_k;
//    std::map<int,float> div3_k;s
    std::vector<float> div_k2;
    for (int i = 1; i < div_k.size(); i++) {
        float k_d2 = (div_k[i - 1] - div_k[i]) / 2;
        div_k2.push_back(k_d2);
//        div3_k.insert(std::make_pair(int(div_k[i]),(k_d2)));
//        div2_k.insert(std::make_pair((k_d2),));
//        k_test << distance_mat_ave[i-1] << "," << k_d2 << std::endl;
    }
    float max22 = -10000;

//
//    std::cout << " key " << (div2_k.begin())->first << " value " << (div2_k.begin())->second << std::endl;
//    std::cout << " key " << (--div2_k.end())->first << " value " << (--div2_k.end())->second << std::endl;
//    std::cout << " key " << (div3_k.begin())->first << " value " << (div3_k.begin())->second << std::endl;
    float knee_point = 0;
}

float DBSCAN::CalcDistance(Point &pt1, Point &pt2) {
    float distance = std::sqrt((pt1.x - pt2.x) * (pt1.x - pt2.x) + (pt1.y - pt2.y) * (pt1.y - pt2.y));
    return distance;
}

//边缘点是否包括进来了？
//corepoint 和 corepoint 是否有join？
bool DBSCAN::ExpandCluster(int pt_idx, int &cluster_idx) {
    std::vector<int> not_core_point;
    std::deque<int> seed_point;
    seed_point.emplace_back(pt_idx);
    for (int col = 0; col < distance_mat_.size(); col++) {
        if (pt_idx == col) { continue; }
        if (distance_mat_[pt_idx][col] <= radius_) {
            seed_point.push_back(col);
        }
    }
    if (seed_point.size() < min_pts_) {
        not_core_point.push_back(pt_idx);
        label_[pt_idx] = -1;
        return false;
    }
    for (int i = 0; i < seed_point.size(); i++) {
        label_[seed_point[i]] = cluster_idx;
    }
    seed_point.pop_front();
    while (!seed_point.empty()) {
        std::vector<int> tmp;
        int row = seed_point.front();
        for (int col = 0; col < distance_mat_.size(); col++) {
            if (distance_mat_[row][col] <= radius_) {
                tmp.push_back(col);
            }
        }
        if (tmp.size() >= min_pts_) {
            for (int i = 0; i < tmp.size(); i++) {
                if (label_[tmp[i]] != 0) {
                    continue;
                }
                seed_point.emplace_back(tmp[i]); // tmp[i]继续作为core点查找
                label_[tmp[i]] = cluster_idx;
            }
        }
        seed_point.pop_front();
    }
    cluster_idx++;
    return true;
}


void DBSCAN::GetCluster() {
    int cluster_idx = 1;
    distance_mat_.resize(all_pt_.size(), std::vector<float>(all_pt_.size(), 0));
    label_ = std::vector<int>(all_pt_.size(), 0);
    //calc distance for each point
    for (int i = 0; i < all_pt_.size(); i++) {
        for (int j = i; j < all_pt_.size(); j++) {
            if (i != j) {
                distance_mat_[i][j] = CalcDistance(all_pt_[i], all_pt_[j]);
                distance_mat_[j][i] = distance_mat_[i][j];
            }
        }
    }
    std::deque<std::deque<int>> core_point(all_pt_.size());
    for (int i = 0; i < all_pt_.size(); i++) {
        if (label_[i] != 0) { continue; }
        ExpandCluster(i, cluster_idx);
    }
    JoinEdgePoint();
    float count = eval();
    std::cout << "[DEBUG] count" << count << std::endl;
}


void DBSCAN::viz() {
    for (int pt = 0; pt < all_pt_.size(); pt++) {
        labels << all_pt_[pt].x << "," << all_pt_[pt].y << "," << label_[pt] << std::endl;
    }

}


class GrahamScan {
public:
    Point s_low_;
    std::vector<Point> all_pt_;

    GrahamScan() = default;

    GrahamScan(std::vector<Point> &all_pt) : all_pt_(all_pt) {};

    ~GrahamScan() {};

    void prioritysort(std::vector<Point> &a);

    void GetPolygen(std::vector<Point> &all_pt, std::vector<Point> &for_esd);

    static float cross2d(const Point &a, const Point &b);

    static float dot2d(const Point &a, const Point &b);

    int to_left(const Point &a, const Point &b, const Point &c);

    void ExtractColinePoint();

    void DoGrahamScan(std::vector<Point> &all_pt, std::vector<Point> &for_esd);

};


void GrahamScan::ExtractColinePoint() {


}


void GrahamScan::prioritysort(std::vector<Point> &all_pt) {
//选取po x最小y最大 （右下脚的目标）
//x 最小和y最大不是同一个点的时候
    for (int i = 0; i < all_pt.size(); i++) {
        std::cout << " sort test1 " << all_pt[i].x << "," << all_pt[i].y << std::endl;
    }
    sort(all_pt.begin(), all_pt.end(), [](const Point &a, const Point &b) -> bool {
        return (a.x<(b.x - 1e-3F) ? true : a.x>(b.x + 1e-3F) ? false : a.y < b.y);
    });
    for (int i = 0; i < all_pt.size(); i++) {
        orig_point2 << i << "," << all_pt[i].x << "," << all_pt[i].y << std::endl;
    }
    s_low_ = all_pt.front();
    auto compare = [&](const Point &a, const Point &b) {
        // OA -OB = BA
        // OA - OC = CA
        // 判断 BA 和 CA的顺序
        Point s_a{a.x - s_low_.x, a.y - s_low_.y};
        Point s_b{b.x - s_low_.x, b.y - s_low_.y};
        float tmp = cross2d(s_a, s_b);
        //如果是基本共线，选取角度小的那个
        if (fabs(tmp) < 1e-6) {
            return dot2d(s_a, s_a) < dot2d(s_b, s_b);
        } else {
            return tmp > 0;
        }
    };
    //以s为基准逆时针排序
    sort(all_pt.begin(), all_pt.end(), compare);
    for (int i = 0; i < all_pt.size(); i++) {
        polygen_test2 << i << "," << all_pt[i].x << "," << all_pt[i].y << std::endl;
    }
}

//连接每一个点，每链接一个点检测连线的走向是否是逆时针，如果是则留下该点的前一个点，反之去除前一个点，使之与前面第二个点直接连接
int GrahamScan::to_left(const Point &p0, const Point &p1, const Point &p2) {
    //OP1 - OP0   = P0P1
    //OP2 - OP1 = P1P2
    //计算P0P1 x P1P2的叉积
    Point line_1{p1.x - p0.x, p1.y - p0.y};
    Point line_2{p2.x - p1.x, p2.y - p1.y};
    float tmp = cross2d(line_1, line_2);
    std::cout << "tmp" << tmp << std::endl;
    //需要再验证下是否=0，是否共线
    if (fabs(tmp) < 1e-6) return 0;
    else if (tmp > 0) {
        return 1;
    } else {
        return 2;
    }
}

void GrahamScan::DoGrahamScan(std::vector<Point> &all_pt, std::vector<Point> &for_esd) {
    prioritysort(all_pt);
    GetPolygen(all_pt, for_esd);
}

void GrahamScan::GetPolygen(std::vector<Point> &all_pt, std::vector<Point> &for_esd) {
    deque<Point> S{all_pt.size()};
    deque<Point> T{all_pt.size()};
    for_esd.reserve(all_pt.size());
    S.resize(0);
    T.resize(0);
    S.push_back(all_pt[0]);
    S.push_back(all_pt[1]);
    for (int i = 2; i < all_pt.size(); i++) {
        T.push_back(all_pt[i]);
    }
    int m = 2;
    for (int i = m; i < all_pt.size(); i++) {
        while (S.size() > 1) {
            if (to_left(*(S.end() - 2), *(S.end() - 1), all_pt[i]) == 1) {
                S.push_back(all_pt[i]);
                break;
            } else {
                S.pop_back();
            }
        }
    }

    for (int i = 0; i < S.size(); i++) {
        final_polygen2 << i << "," << S[i].x << "," << S[i].y << std::endl;
    }
}


float GrahamScan::cross2d(const Point &a, const Point &b) {
    return a.x * b.y - b.x * a.y;
}

float GrahamScan::dot2d(const Point &a, const Point &b) {
    return a.x * b.x + a.y * b.y;
}

//以空格为分隔符，逐行读取并划分
std::vector<Point> ReadTxt2(string &file) {
    string line;
    std::vector<Point> pp;
    ifstream infile(file);
    while (getline(infile, line, '\n')) {
        stringstream ss(line);
        Point p;
        ss >> p.x >> p.y;
        pp.push_back(p);
    }
    return pp;
}

std::vector<float> ReadTxt(string &file) {
    string line;
    ifstream infile(file);
    std::vector<string> x;

    std::vector<float> x1;
    if (infile) {
        while (getline(infile, line, '\n')) {
            x.push_back(line);
        }
        for (int i = 0; i < x.size(); i++) {
            x1.push_back(stof(x[i]));
        }
        return x1;
    }

    return x1;
}

int main() {
    float eps = 0.25;
    int min_pts = 5;
    float eps2 = 30;
    int min_pts2 = 5;
    std::vector<float> x1;
    std::vector<float> x2;
    std::vector<float> y1;
    std::vector<float> y2;
    std::vector<float> cx1;
    std::vector<float> cy1;
    std::vector<Point> ped;
    string file_x1 = "/Users/xinyi.wang2/CLionProjects/untitled/x1.txt";
    string file_x2 = "/Users/xinyi.wang2/CLionProjects/untitled/x2.txt";
    string file_y1 = "/Users/xinyi.wang2/CLionProjects/untitled/y1.txt";
    string file_y2 = "/Users/xinyi.wang2/CLionProjects/untitled/y2.txt";
    string file_cx1 = "/Users/xinyi.wang2/CLionProjects/untitled/cx1.txt";
    string file_cy1 = "/Users/xinyi.wang2/CLionProjects/untitled/cy1.txt";
    string file_ped = "/Users/xinyi.wang2/CLionProjects/untitled/ped.txt";
    x1 = ReadTxt(file_x1);
    x2 = ReadTxt(file_x2);
    y1 = ReadTxt(file_y1);
    y2 = ReadTxt(file_y2);
    cx1 = ReadTxt(file_cx1);
    cy1 = ReadTxt(file_cy1);
    ped = ReadTxt2(file_ped);


    DBSCAN deb3(ped, 3);
    deb3.GetCluster();
    deb3.viz();
    std::vector<Point> c_all_pts;
    std::vector<Point> all_pts;
    std::cout << " c_all_pts " << c_all_pts.size() << std::endl;
    for (int i = 0; i < cx1.size(); i++) {
        Point pt;
        pt.x = cx1[i];
        pt.y = cy1[i];
        c_all_pts.emplace_back(pt);
    }

//    std::vector<Point> test_pt2;
//    test_pt2.resize(30);
//    std::copy( c_all_pts.begin(), c_all_pts.begin() + 10, test_pt2.begin());




    DBSCAN deb2(c_all_pts, min_pts2);
    deb2.GetCluster();
    deb2.viz();

    std::vector<Point> for_esd;
    GrahamScan grahamScan2;
    grahamScan2.prioritysort(c_all_pts);
    grahamScan2.GetPolygen(c_all_pts, for_esd);
    grahamScan2.ExtractColinePoint();
    for (int i = 0; i < x1.size(); i++) {
        Point pt;
        pt.x = x1[i];
        pt.y = y1[i];
        all_pts.emplace_back(pt);
        Point pt2;
        pt2.x = x2[i];
        pt2.y = y2[i];
        all_pts.emplace_back(pt2);
    }

    std::vector<Point> test_pt;
    test_pt.resize(30);
    std::copy(all_pts.begin(), all_pts.begin() + 30, test_pt.begin());
    std::vector<Point> for_esd2;
    DBSCAN deb(all_pts, min_pts, eps);
    deb.GetCluster();
    deb.viz();

    DBSCAN dbtest(all_pts, min_pts);
    dbtest.GetCluster();
    dbtest.viz();

    GrahamScan grahamScan;
    grahamScan.prioritysort(test_pt);
    grahamScan.GetPolygen(test_pt, for_esd);


    return 0;
}
