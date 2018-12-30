crop: 300 180 300 300


./hw1 -input img/yaoming.jpg -quentize 1 -output quantize_1.png
./hw1 -input img/yaoming.jpg -FloydSteinbergDither 1 -output Fdither_1.png


//
0 -> IMAGE_SAMPLING_POINT,
1 -> IMAGE_SAMPLING_BILINEAR,
2 ->IMAGE_SAMPLING_GAUSSIAN
//
./hw1 -input img/tiger.bmp -sampling 0 -scale 1.2 1.2 -output scalep.png
./hw1 -input img/tiger.bmp -sampling 1 -scale 1.2 1.2 -output scaleb.png
./hw1 -input img/tiger.bmp -sampling 2 -scale 1.2 1.2 -output scaleg.png

./hw1 -input img/cube.jpg -sampling 0 -rotate 90 -output rotatep.png
./hw1 -input img/jacky.jpg -sampling 1 -rotate 90 -output rotateb.png
./hw1 -input img/jacky.jpg -sampling 2 -rotate 90 -output rotateg.png
