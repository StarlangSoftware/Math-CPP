Video Lectures
============

[<img src=https://github.com/StarlangSoftware/Math/blob/master/video.jpg width="50%">](https://youtu.be/GhcoaVi0SMs)

For Developers
============
You can also see [Java](https://github.com/starlangsoftware/Math), [Python](https://github.com/starlangsoftware/Math-Py), [Cython](https://github.com/starlangsoftware/Math-Cy), [Swift](https://github.com/starlangsoftware/Math-Swift), [C](https://github.com/starlangsoftware/Math-C), [Js](https://github.com/starlangsoftware/Math-Js), or [C#](https://github.com/starlangsoftware/Math-CS) repository.

## Requirements

* [C++ Compiler](#cpp)
* [Git](#git)


### CPP
To check if you have compatible C++ Compiler installed,
* Open CLion IDE 
* Preferences >Build,Execution,Deployment > Toolchain  

### Git

Install the [latest version of Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

## Download Code

In order to work on code, create a fork from GitHub page. 
Use Git for cloning the code to your local or below line for Ubuntu:

	git clone <your-fork-git-link>

A directory called Math-CPP will be created. Or you can use below link for exploring the code:

	git clone https://github.com/starlangsoftware/Math-CPP.git

## Open project with CLion IDE

To import projects from Git with version control:

* Open CLion IDE , select Get From Version Control.

* In the Import window, click URL tab and paste github URL.

* Click open as Project.

Result: The imported project is listed in the Project Explorer view and files are loaded.


## Compile

**From IDE**

After being done with the downloading and opening project, select **Build Project** option from **Build** menu. After compilation process, user can run Math-CPP.

Detailed Description
============

+ [Vector](#vector)
+ [Matrix](#matrix)
+ [Distribution](#distribution)

## Vector

Bir vektör yaratmak için:

	a = Vector(ArrayList<Double> values)

Vektörler eklemek için

	void add(Vector v)

Çıkarmak için

	void subtract(Vector v)
	Vector difference(Vector v)

İç çarpım için

	double dotProduct(Vector v)
	double dotProduct()

Bir vektörle cosinüs benzerliğini hesaplamak için

	double cosineSimilarity(Vector v)

Bir vektörle eleman eleman çarpmak için

	Vector elementProduct(Vector v)

## Matrix

3'e 4'lük bir matris yaratmak için

	a = Matrix(3, 4)

Elemanları rasgele değerler alan bir matris yaratmak için

	Matrix(int row, int col, double min, double max)

Örneğin, 

	a = Matrix(3, 4, 1, 5)
 
3'e 4'lük elemanları 1 ve 5 arasında değerler alan bir matris yaratır.

Birim matris yaratmak için

	Matrix(int size)

Örneğin,

	a = Matrix(4)

4'e 4'lük köşegeni 1, diğer elemanları 0 olan bir matris yaratır.

Matrisin i. satır, j. sütun elemanını getirmek için 

	double getValue(int rowNo, int colNo)

Örneğin,

	a.getValue(3, 4)

3. satır, 4. sütundaki değeri getirir.

Matrisin i. satır, j. sütunundaki elemanı değiştirmek için

	void setValue(int rowNo, int colNo, double value)

Örneğin,

	a.setValue(3, 4, 5)

3. satır, 4.sütundaki elemanın değerini 5 yapar.

Matrisleri toplamak için

	void add(Matrix m)

Çıkarmak için 

	void subtract(Matrix m)

Çarpmak için 

	Matrix multiply(Matrix m)

Elaman eleman matrisleri çarpmak için

	Matrix elementProduct(Matrix m)

Matrisin transpozunu almak için

	Matrix transpose()

Matrisin simetrik olup olmadığı belirlemek için

	boolean isSymmetric()

Determinantını almak için

	double determinant()

Tersini almak için

	void inverse()

Matrisin eigenvektör ve eigendeğerlerini bulmak için

	ArrayList<Eigenvector> characteristics()

Bu metodla bulunan eigenvektörler eigendeğerlerine göre büyükten küçüğe doğru 
sıralı olarak döndürülür.

## Distribution

Verilen bir değerin normal dağılımdaki olasılığını döndürmek için

	static double zNormal(double z)

Verilen bir olasılığın normal dağılımdaki değerini döndürmek için

	static double zInverse(double p)

Verilen bir değerin chi kare dağılımdaki olasılığını döndürmek için

	static double chiSquare(double x, int freedom)

Verilen bir olasılığın chi kare dağılımdaki değerini döndürmek için

	static double chiSquareInverse(double p, int freedom)

Verilen bir değerin F dağılımdaki olasılığını döndürmek için

	static double fDistribution(double F, int freedom1, int freedom2)

Verilen bir olasılığın F dağılımdaki değerini döndürmek için

	static double fDistributionInverse(double p, int freedom1, int freedom2)

Verilen bir değerin t dağılımdaki olasılığını döndürmek için

	static double tDistribution(double T, int freedom)

Verilen bir olasılığın t dağılımdaki değerini döndürmek için

	static double tDistributionInverse(double p, int freedom)
