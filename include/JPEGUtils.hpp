#include <vector>

class JPEGImage {
 public:
  JPEGImage();
  JPEGImage(std::vector<std::pair<int, int>> sizes);
  ~JPEGImage();

  std::pair<int, int> GetComponentShape(int component);

  int& at(int row, int col, int component);

 private:
  std::vector<int*> image_components_;
  std::vector<std::pair<int, int>> components_shape;
};