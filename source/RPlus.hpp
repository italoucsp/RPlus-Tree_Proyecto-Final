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

      Entry(const std::shared_ptr<RData_type> record) {
        record_ = record;
        mbr_ = record();
      }

      Entry(const std::shared_ptr<RPNode> son_ptr) {
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

      std::shared_ptr<RPNode> get_son() const noexcept{
        return son_ptr_;
      }

    private:
      KDRect<K_Dimensions> mbr_;
      std::shared_ptr<RPNode> son_ptr_;
      std::shared_ptr<RData_type> record_;
    };

    struct RPNode {
      typedef Entry* iterator;
      typedef const Entry* const_iterator;

      iterator begin() { return fields.begin(); }

      iterator end() { return last_; }

      const_iterator begin() const { return fields.begin(); }

      const_iterator end() const { return last_; }

      RPNode(std::size_t level = 0) { 
        last_ = fields.before_begin();
        size_ = 0;
        level_ = level;
      }

      std::size_t size() const noexcept{
        return size_;
      }

      std::size_t get_level() const noexcept {
        return level_;
      }

      bool is_leaf() noexcept { return level_ == 0; }

      bool is_overflowed() noexcept { return size_ > Node_Size; }

      void insert(const Entry entry) {
        last_ = fields.insert_after(last_, entry);
        ++size_;
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

      std::shared_ptr<RPNode> split(std::size_t& axis, double& cutline) {
        std::shared_ptr<RPNode> other_half = std::make_shared<RPNode>(level_);
        std::forward_list<Entry> left_group, right_group;
        for (Entry& entry : fields) {
          entry.get_mbr();
        }
        return other_half;
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
      std::shared_ptr<RPNode> cnode = choose_leaf(data(), ancestors);
      std::shared_ptr<RPNode> splitted_node_left, splitted_node_right;
      cnode->insert(Entry(std::make_shared<RData_type>(data)));
      ancestors.push(cnode);
      while (ancestors.top()->is_overflowed()) {
        std::size_t current_axis;
        double current_cutline;
        ancestors.top()->find_best_partition(current_axis, current_cutline);
        splitted_node_right = ancestors.top()->split(current_axis, current_cutline);
        splitted_node_left = ancestors.top();
        ancestors.pop();
        if (!ancestors.empty()) {//normal split-insertion operation
          ancestors.top()->insert(splitted_node_left);
          ancestors.top()->insert(splitted_node_right);
        }
        else {//new root by split-insertion operation
          std::shared_ptr<RPNode> newroot = std::make_shared<RPNode>(root_->get_level());
          newroot->insert(splitted_node_left);
          newroot->insert(splitted_node_right);
          root_ = newroot;
          return;
        }
      }
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
    std::shared_ptr<RPNode> choose_leaf(const RContainer_type& val_container,
                                std::stack<std::shared_ptr<RPNode>>& ancestors_path) {
      std::shared_ptr<RPNode> cnode = root;
      while (!cnode->is_leaf()) {
        ancestors_path.push(cnode);
        std::shared_ptr<RPNode> temp = cnode;
        for (Entry& entry : temp) {
          cnode = entry.get_son();
          if (entry.get_mbr().overlaps(val_container))
            break;
        }
      }
      return cnode;
    }
  };
}

#endif //R_PLUS_HPP