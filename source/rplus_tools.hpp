#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
//A-Z

//This file only contains tools for R+

void console_output__LogTime() {

}

template<size_t K_Dimensions>
class KDGeometry {
protected:
  static std::size_t KDGeometry_type() { return typeid(KDGeometry<K_Dimensions>).hash_code(); }
public:
  static const std::size_t Dimensionality = K_Dimensions;
  template<typename Container, std::size_t Dimensionality>
  friend class KDRecord;
};

template<size_t K_Dimensions>
class KDPoint : public KDGeometry<K_Dimensions> {
  std::array<double, K_Dimensions> axis_values_;
public:
  KDPoint() {
    axis_values_.fill(double(0));
  }

  double operator[](size_t idx) const {
    return axis_values_[idx];
  }

  KDPoint<K_Dimensions>& operator=(const KDPoint<K_Dimensions>& point_value) {
    std::size_t idx = 0;
    for (double& value : axis_values_) {
      value = elem_value.axis_values_[idx++];
    }
  }

  static KDPoint get_max() {
    KDPoint maxpoint;
    for (double& value : maxpoint.axis_values_) {
      value = std::numeric_limits<double>::max();
    }
    return maxpoint;
  }

  static KDPoint get_min() {
    KDPoint minpoint;
    for (double& value : minpoint.axis_values_) {
      value = std::numeric_limits<double>::min();
    }
    return minpoint;
  }

  template<size_t K_Dimensions, typename stream_input_type>
  friend stream_input_type& operator>>(stream_input_type& is, KDPoint<K_Dimensions>& point);
};

template<size_t K_Dimensions, typename stream_input_type>
stream_input_type& operator>>(stream_input_type& is, KDPoint<K_Dimensions>& point) {
  for (double& value : point.axis_values_)
    is >> value;
  return is;
}

template<size_t K_Dimensions>
class KDRect : public KDGeometry<K_Dimensions> {
  KDPoint<K_Dimensions> bottom_left_, top_right_;
public:
  KDRect() {
    bottom_left_ = KDPoint<K_Dimensions>::get_max();
    top_right_ = KDPoint<K_Dimensions>::get_min();
  }

  KDPoint<K_Dimensions> get_bl() { return bottom_left; }

  KDPoint<K_Dimensions> get_tr() { return top_right; }

  KDRect<K_Dimensions>& operator=(const KDPoint<K_Dimensions>& point_value) {
    bottom_left_ = point_value;
    top_right_ = point_value;
  }

  void enlarge(const KDRect<K_Dimensions>& other) {
    for (size_t idx(0); idx < K_Dimensions; ++idx) {
      bottom_left_[idx] = std::min(bottom_left_[idx], other[idx]);
      top_right_[idx] = std::max(top_right_[idx], other[idx]);
    }
  }

  bool overlaps(const KDRect<K_Dimensions>& rect) {
    for (size_t idx(0); idx < K_Dimensions; ++idx) {
      if (bottom_left_[idx] > rect.top_right_[idx] ||
        top_right_[idx] < rect.bottom_left_[idx])
        return false;
    }
    return true;
  }

  bool overlaps(const KDPoint<K_Dimensions>& point) {
    for (size_t idx(0); idx < K_Dimensions; ++idx) {
      if (bottom_left_[idx] > point[idx] || top_right_[idx] < point[idx])
        return false;
    }
    return true;
  }
};

template<typename Container, std::size_t Dimensionality = Container::Dimensionality>
class KDRecord {
  typedef KDPoint<Dimensionality> Container_t_Point;
  typedef KDRect<Dimensionality> Container_t_Rect;
public:
  static const std::size_t RDimensionality = Dimensionality;
  typedef Container RContainer;

  virtual RContainer operator()() = 0;
  static bool check_container_class() {
    return (Container::KDGeometry_type() == Container_t_Point::KDGeometry_type() ||
      Container::KDGeometry_type() == Container_t_Rect::KDGeometry_type());
  }
};
/*
using namespace std;

#define ALERT(message) cerr << "[ERROR] : " << message << endl;
#define SAY(message) cout << "[STEP] : " << message << endl;
#define ERROR_M_N_VALUES "The number of dimensions and/or the number of entries per node should be greater than 1."
#define ERROR_FF_VALUE "The value for fill factor should be between 2 and M."
#define ERROR_NODE_OFR "The index is out of range in the node."
#define ERROR_EMPTY_TREE "This R+ Tree is empty."

const size_t T_DIMENSIONS_NUM = 19;//Number of dimensions in the dataset
const size_t KUSED_DIMENSIONS = 14;//Number of dimensions that will be used

const char csv_delimiter = ';';

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HyperPoint : DATA or Bound for HyperRectangle

template<typename T, size_t N>
struct HyperPoint {
  HyperPoint();
  HyperPoint(array<T, N> data);
  HyperPoint(array<T, N> data, string sg_name);
  HyperPoint<T, N>& operator=(const HyperPoint<T, N> &other);
  T& operator[](size_t index);
  T operator[](size_t index) const;
  string get_songs_name();
  void show_data();

private:
  array<T, N> multidata;
  string songs_name;
};

template<typename T, size_t N>
HyperPoint<T, N>::HyperPoint() {
  multidata.fill(T(0));
}

template<typename T, size_t N>
HyperPoint<T, N>::HyperPoint(array<T, N> data) {
  for (size_t i(0); i < N; ++i) {
    multidata[i] = data[i];
  }
}

template<typename T, size_t N>
HyperPoint<T, N>::HyperPoint(array<T, N> data, string sg_name) {
  for (size_t i(0); i < N; ++i) {
    multidata[i] = data[i];
  }
  songs_name = sg_name;
}

template<typename T, size_t N>
HyperPoint<T, N>& HyperPoint<T, N>::operator=(const HyperPoint<T, N>& other) {
  songs_name = other.songs_name;
  for (size_t i(0); i < N; ++i)
    multidata[i] = other.multidata[i];
  return *this;
}

template<typename T, size_t N>
T& HyperPoint<T, N>::operator[](size_t index) {
  return multidata[index];
}

template<typename T, size_t N>
T HyperPoint<T, N>::operator[](size_t index) const{
  return multidata[index];
}

template<typename T, size_t N>
string HyperPoint<T, N>::get_songs_name() {
  return songs_name;
}

template<typename T, size_t N>
void HyperPoint<T, N>::show_data() {
  cout << "\tHyperPoint<" << N << "> : {"; for (size_t i(0); i < N; ++i) cout << setprecision(10) << multidata[i] << ((i != N - 1) ? "," : ""); cout << "}\n";
}

template<typename T, size_t N>
istream& operator>>(istream& is, HyperPoint<T, N>& hpoint) {
  array<T, N> raw_data;
  for(size_t i(0); i < N; ++i)
    is >> raw_data[i];
  hpoint = HyperPoint<T, N>(raw_data);
  return is;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Class for the MBR of each entry/node in the R+
template<typename T, size_t N>
struct HyperRectangle {
  HyperRectangle();
  HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B);
  HyperRectangle<T, N>& operator=(const HyperRectangle<T, N>& other);
  bool overlaps(const HyperRectangle<T, N> &other);
  bool contains(const HyperPoint<T, N> &point);
  void adjust(const HyperRectangle<T, N> &other);
  pair<HyperPoint<T, N>, HyperPoint<T, N>> get_boundaries();
  double get_hypervolume();
  void show_rect();

  template<typename T, size_t N>
  friend HyperRectangle<T, N> make_hyper_rect(HyperPoint<T, N> &h_point);

private:
  HyperPoint<T, N> bottom_left, top_right;
};

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle() {
  //...
}

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B) {
  for (size_t i(0); i < N; ++i) {
    bottom_left[i] = min(A[i], B[i]);
    top_right[i] = max(A[i], B[i]);
  }
}

template<typename T, size_t N>
HyperRectangle<T, N>& HyperRectangle<T, N>::operator=(const HyperRectangle<T, N>& other) {
  bottom_left = other.bottom_left;
  top_right = other.top_right;
  return *this;
}

template<typename T, size_t N>
bool HyperRectangle<T, N>::overlaps(const HyperRectangle<T, N> &other) {
  for (size_t i(0); i < N; ++i) {
    if (!(bottom_left[i] <= other.top_right[i] &&
          top_right[i] >= other.bottom_left[i])) {
      return false;
    }
  }
  return true;
}

template<typename T, size_t N>
bool HyperRectangle<T, N>::contains(const HyperPoint<T, N> &point) {
  for (size_t i(0); i < N; ++i) {
    if (!(bottom_left[i] <= point[i] && point[i] <= top_right[i]))
      return false;
  }
  return true;
}

template<typename T, size_t N>
void HyperRectangle<T, N>::adjust(const HyperRectangle<T, N> &other) {
  for (size_t i(0); i < N; ++i) {
    bottom_left[i] = min(other.bottom_left[i], bottom_left[i]);
    top_right[i] = max(other.top_right[i], top_right[i]);
  }
}

template<typename T, size_t N>
pair<HyperPoint<T, N>, HyperPoint<T, N>> HyperRectangle<T, N>::get_boundaries() {
  return make_pair(bottom_left, top_right);
}

template<typename T, size_t N>
double HyperRectangle<T, N>::get_hypervolume() {
  double hypervolume = 1.0;
  for (size_t i(0); i < N; ++i) {
    hypervolume *= double(top_right[i] - bottom_left[i]);
  }
  return hypervolume;
}

template<typename T, size_t N>
void HyperRectangle<T, N>::show_rect() {
  cout << "\t\tHyperRectangle<" << N << "> :   bottom_left: "; bottom_left.show_data();
  cout << "\t\t\t\t\ttop_right: "; top_right.show_data();
}

template<typename T, size_t N>
HyperRectangle<T, N> make_hyper_rect(HyperPoint<T, N> &h_point) {
  HyperRectangle<T, N> conv;
  conv.bottom_left = conv.top_right = h_point;
  return conv;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//CSV file reader : path of the file | features that were considered | id(name of the song) | container for the data in hypepoints
template<typename T, size_t N>
void read_data_from_file(string file_path, vector<string> &features, string id, vector<HyperPoint<T, N>> &db_container) {
  string cols_labels, row_data_line;
  ifstream data_set_file(file_path);
  if (!data_set_file.is_open()) {
    ALERT("Couldn\'t open the file " + file_path)
    exit(1);
  }
  if (N + 1 < features.size()) {
    ALERT("The value of N is not enough for the given features.")
    exit(1);
  }

  getline(data_set_file, cols_labels);
  string label;
  istringstream iss_cols_labels(cols_labels);
  array<bool, T_DIMENSIONS_NUM> table_check_features;
  table_check_features.fill(false);
  size_t tcfi(0), id_id(0);
  while (!iss_cols_labels.eof()) {
    getline(iss_cols_labels, label, csv_delimiter);
    for (size_t fi(0); fi < features.size(); ++fi) {
      if (features[fi] == label) {
        table_check_features[tcfi] = true;
        if (label == id)
          id_id = tcfi;
      }
    }
    ++tcfi;
  }
  while (getline(data_set_file, row_data_line)) {
    istringstream iss_cols_data(row_data_line);
    size_t i(0), ii(0);
    string data_col_row, songs_name;
    array<T, KUSED_DIMENSIONS> raw_multidimensional_data;
    while (!iss_cols_data.eof()) {
      getline(iss_cols_data, data_col_row, csv_delimiter);
      if (i < T_DIMENSIONS_NUM && table_check_features[i]) {
        if (id_id == i)
          songs_name = data_col_row;
        else {
          stringstream ss_raw_datacol(data_col_row);
          ss_raw_datacol >> raw_multidimensional_data[ii];
          ++ii;
        }
      }
      ++i;
    }
    HyperPoint<T, KUSED_DIMENSIONS> multidimensional_data(raw_multidimensional_data, songs_name);
    db_container.push_back(multidimensional_data);
  }

  data_set_file.close();
}*/