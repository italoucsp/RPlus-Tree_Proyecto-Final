#include <RPlusTree.hpp>
#include "RPlus.hpp"

/*
-----------------------[ADVANCED DATA STRUCTURES]------------------------
-----------------SPOTIFY DATASET 1921-2020, 160k+TRACKS------------------
-------------------------DATA-STRUCT.: R+ Tree---------------------------
------------------------Italo Mamani Huaricallo--------------------------
---------------------------Arequipa-Perú-2020----------------------------
*/

struct SpotifySongData : public KDRecord<KDPoint<14>> {
  std::string song_name_, artist_, id_;
  double acousticness_, danceability_, energy_, instrumentalness_,
         liveness_, loudness_, speechiness_, tempo_, valence_;
  int64_t duration_ms_;
  bool explicit_, mode_;
  uint16_t key_;
  int16_t year_, release_date;
  uint8_t popularity_;

  KDPoint<14> operator()() {
    KDPoint<14> kdGeoContainer;
    std::stringstream data_setter;
    data_setter << acousticness_ << ' ' << danceability_ << ' ' << 
      duration_ms_ << ' ' << energy_ << ' ' << explicit_ << ' ' << 
      instrumentalness_ << ' ' << key_ << ' ' << liveness_ << ' ' << 
      loudness_ << ' ' << mode_ << ' ' << popularity_ << ' ' << 
      speechiness_ << ' ' << tempo_ << ' ' << valence_;
    data_setter >> kdGeoContainer;
    return kdGeoContainer;
  }
};

int main() {

  ads::RPlusTree<5, 4, SpotifySongData> demo;
  return 0;
}