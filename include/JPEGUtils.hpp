#include <vector>

class JPEGImage {
 public:
  JPEGImage(std::vector<std::pair<int, int>> sizes);
  ~JPEGImage();

  std::pair<int, int> GetComponentShape(int component);

  int& at(int col, int row, int component);

 private:
  std::vector<int*> image_components_;
  std::vector<std::pair<int, int>> components_shape;
}