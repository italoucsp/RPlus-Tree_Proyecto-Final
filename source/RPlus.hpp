#ifndef SOURCE_R_PLUS_HPP
#define SOURCE_R_PLUS_HPP

#include <assert.h>
#include <chrono>
#include <forward_list>
#include <functional>
#include <memory>
#include <typeinfo>
#include <stack>
#include <stdexcept>
#include <stdlib.h>
#ifdef _WIN32
#include <Windows.h>
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#define set_text_color(Text)                                                                                                                      \
  SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), Text);
#else
#define set_text_color(Text) (void(0))
#endif


#include "rplus_tools.hpp"

enum console_colors { COLOR_ERROR = 12, COLOR_WARNING = 14, COLOR_NORMAL = 15 };

#ifdef _DEBUG
      #define Assert_expression(expression, typerr,message) assert(expression)
#else
      #define Assert_expression(expression, typerr, message)                                                                                      \
              try { if (!(expression)) throw typerr(message); }                                                                                   \
              catch (const std::exception& error_) {                                                                                              \
                set_text_color(COLOR_ERROR);                                                                                                      \
                std::cerr << "[ERROR]: " << error_.what() << '\n';                                                                                \
                set_text_color(COLOR_NORMAL);                                                                                                     \
                exit(EXIT_FAILURE);                                                                                                               \
              }
#endif

#define Throw_warning_if_not(expression, message)                                                                                                 \
              try { if (!(expression)) throw std::string(message); }                                                                              \
              catch (const std::string adv_) {                                                                                                    \
                set_text_color(COLOR_WARNING);                                                                                                    \
                std::cerr << "[WARNING]: " << adv_ << '\n';                                                                                       \
                set_text_color(COLOR_NORMAL);                                                                                                     \
              }

#define POINT_APPROACH 1
#define VOLUMETRIC_APPROACH 0

namespace ads {
  namespace {
    typedef std::invalid_argument err_iar;
    typedef std::out_of_range     err_oor;
    typedef std::logic_error      err_log;
  }

  template<std::size_t Node_Size, std::size_t Fill_Factor, typename RData_type, std::size_t K_Dimensions = RData_type::RDimensionality>
  class RPlusTree {
    typedef long double Cost_type;
    typedef typename RData_type::RContainer RContainer_type;
    struct RPNode;
    struct Entry;

    std::function<Cost_type(std::vector<Entry>&, std::size_t, double)> cost_Function_;

    struct Entry {
      Entry() {}

      Entry(const std::shared_ptr<RData_type>& record) {
        record_ = record;
        mbr_ = record();
      }

      Entry(const std::shared_ptr<RPNode>& son_ptr) {
        son_ptr_ = son_ptr;
        mbr_ = son_ptr->calculate_mbr();
      }

      Entry(const Entry& other) {
        mbr_ = other.mbr_;
        son_ptr_ = other.son_ptr_;
        record_ = other.record_;
      }

      KDRect<K_Dimensions>& get_mbr() const noexcept {
        return mbr_;
      }
    private:
      KDRect<K_Dimensions> mbr_;
      std::shared_ptr<RPNode> son_ptr_;
      std::shared_ptr<RData_type> record_;
    };
    struct RPNode {
      RPNode() { last_ = fields.before_begin(); }

      std::size_t size() const noexcept{
        return size_;
      }

      bool is_leaf() { return level_ == 0; }

      void insert(const Entry& entry) {
        last_ = fields.insert_after(last_, entry);
        ++size_;
      }
      
      const Entry& scan_node(const RContainer_type& container_value) {
        for (Entry& entry : fields) {
          if (entry.get_mbr().overlaps(container_value))
            return entry;
        }
      }

      const KDRect<K_Dimensions> calculate_mbr() {
        KDRect<K_Dimensions> ans;
        for (Entry& field : fields) {
          ans.enlarge(field.get_mbr());
        }
        return ans;
      }

      Cost_type sweep(std::vector<Entry>& test_group, std::size_t axis, double& cutline) {
        std::partial_sort(test_group.begin(), test_group.begin() + Fill_Factor, test_group.end(),
          [=](Entry one_, Entry another_) {
          return one_.get_mbr().get_bl()[axis] < another_.get_mbr().get_bl()[axis]; 
        });
        std::vector<Entry> ff_sorted(test_group.begin(), test_group.begin() + Fill_Factor);
        cutline = (*ff_sorted.end()).get_mbr().get_tr()[axis];
        return cost_Function_(ff_sorted, axis, cutline);
      }

      void find_best_partition(std::size_t& axis, double& cutline) {
        Cost_type min_cost = std::numeric_limits<Cost_type>::max();
        std::vector<Entry> to_test(fields.begin(), fields.end());
        for (std::size_t axis_idx(0); axis_idx < K_Dimensions; ++axis_idx) {
          double new_cutline;
          Cost_type new_cost = sweep(to_test, axis_idx, new_cutline);
          if (new_cost < min_cost) {
            axis = axis_idx;
            cutline = new_cutline;
            min_cost = new_cost;
          }
        }
      }
    private:
      std::forward_list<Entry> fields;
      typename std::forward_list<Entry>::iterator last_;
      std::size_t size_, level_;
    };
    
    std::shared_ptr<RPNode> root_;
  public://public methods
    explicit RPlusTree() {
      Assert_expression(RData_type::check_container_class(), err_iar,
        "The given type for container class can not be used, only KDRect or KDPoint.");
      Assert_expression(RData_type::RDimensionality < 20, err_oor,
        "The given number of dimensions value is too long.");
      Assert_expression(RData_type::RDimensionality > 1, err_oor,
        "The given number of dimensions value should be greater than 1.");
      Assert_expression(Fill_Factor < Node_Size, err_log,
        "The given value for fill factor should be less than node size's value");
      root_ = std::make_shared<RPNode>();
      std::cout << "Arbol creado" << std::endl;
    }
    virtual ~RPlusTree() {
      root_.reset();
    }
    void insert(const RData_type& data) {
      std::stack<std::shared_ptr<RPNode>> ancestors;
      std::shared_ptr<RPNode> cnode = choose_leaf(ancestors);
    }
    void assign(const std::vector<RData_type>& data_set) {
      std::chrono::time_point<std::chrono::high_resolution_clock> start_time, end_time;
      start_time = std::chrono::high_resolution_clock::now();
      for (RData_type& data : data_set) {
        insert(data);
      }
      end_time = std::chrono::high_resolution_clock::now();
    }
    std::vector<RData_type> knn_query(std::size_t k, KDPoint<K_Dimensions> center) {
      std::chrono::time_point<std::chrono::high_resolution_clock> start_time, end_time;
      start_time = std::chrono::high_resolution_clock::now();
      //query
      end_time = std::chrono::high_resolution_clock::now();

    }
  private://private methods
    std::shared_ptr<RPNode> choose_leaf(std::stack<std::shared_ptr<RPNode>>& ancestors_path) {
      std::shared_ptr<RPNode> cnode = root;
      return nullptr;
    }
  };
}

#endif //R_PLUS_HPP