# C++用線形代数ライブラリ「libLinalg」

## 全体的に準備中

## 演算について
普通に``+``，``-``，``*``，``/``で，ベクトル同士や，行列同士，スカラー倍などの演算が可能

## 関数について

### Vector classのメンバ関数
- ``void print()``<br>ベクトルを表示する関数
- ``setLength(const int length)``<br>ベクトルの長さを``length``に設定する関数
- ``int length()``<br>ベクトルの長さを返す関数
- ``Vector<U> cast<U>()``<br>ベクトルをcastする関数
- ``substitute(const int index, const T number_to_substitute)``<br>ベクトルの``index``番目に``number_to_substitute``を代入する関数
- ``substitute(std::vector<T> vector)``<br>ベクトルに``vector``を代入する関数
- ``int at(const int index)``<br>ベクトルの``index``番目の要素を取ってくる関数
- ``bool isZero()``<br>ベクトルが零ベクトルであるかを判定する関数
- ``void setSeed(const unsigned int seed)``<br>乱数のシードを設定する関数
- ``void random()``<br>ベクトルをランダムに生成する関数
- ``void random(const double lower_bound, const double upper_bound)``<br>ベクトルをランダムに生成する関数．但し，各要素は``lower_bound``以上``upper_bound``以下である
- ``T innerProduct(Vector<T> vector)``<br>``vector``との標準内積を求める関数
- ``norm(const unsigned int p = 2)``<br>$l_p$-ノルムを求める関数

### Matrix classのメンバ関数
- ``Matrix<T> power(const int n)``<br>行列の``n``乗を計算する関数
- ``void print()``<br>行列を表示する関数
- ``void print(const char *print_mode)``<br>行列を表示する関数で，表示の方法を指定できる．今は``python``，``c``，``TeX``がある
- ``T at(const int row_index, const int column_index)``<br>行列の(``row_index``, ``column_index``)-成分を返す関数
- ``void setSeed(const unsigned int seed)``<br>乱数のシードを設定する関数
- ``void random()``<br>行列をランダムに生成する関数
- ``void random(const double lower_bound, const double upper_bound)``<br>行列をランダムに生成する関数．但し，各要素は``lower_bound``以上``upper_bound``以下である
- ``void setDimensions(const int num_rows, num int num_columns)``<br>行列のサイズを``num_rows``×``num_columns``に設定する関数
- ``void setIdentity(const int size)``<br>行列を``size``次単位行列に設定する関数
- ``void setDiagonal(Vector<T> diagonal_vector)``<br>対角成分が``diagonal_vector``であるような対角行列に設定する関数
- ``bool isRegular()``<br>行列が正則かどうか判定する関数
- ``bool isSquare()``<br>行列が正方かどうか判定する関数
- ``bool isDiagonal()``<br>行列が対角行列かどうか判定する関数
- ``bool isUniModular()``<br>行列がユニモジュラかどうか判定する関数
- ``void substitute(const std::vector<std::vector<T>> matrix)``<br>行列に``matrix``を代入する関数
- ``void substituteRow(const int index, std::vector<T> vector)``
- ``void substituteRow(const int index, Vector<T> vector)``<br>行列の``index``行目にベクトル``vector``を代入する関数
- ``void substituteColumn(const int index, std::vector<T> vector)``
- ``void substituteColumn(const int index, Vector<T> vector)``<br>行列の``index``列目にベクトル``vector``を代入する関数
- ``void substitute(const int row_index, const int column_index, const T number)``<br>行列の(``row_index``，``column_index``)-成分に``number``を代入する関数
- ``Matrix<U> cast<U>()``<br>行列をcastする関数
- ``Vector<T> row(const int index)``<br>行列の第``index``行を取り出す関数
- ``int rows()`<br>行列の行数を返す関数
- ``int columns()``<br>行列の列数を返す関数
- ``T det()``
- ``T determinant()``<br>行列の行列式を求める関数
- ``Matrix<T> transpose()``
- ``Matrix<T> trans()``<br>転置行列を求める関数
- ``T trace()``<br>行列のトレースを求める関数
- ``T coFactor(const int row_index, const int column_index)``<br>行列の(``row_index``, ``column_index``)余因子を求める関数
- ``Matrix<T> adjugate()``<br>行列の余因子行列を求める関数
- ``Matrix<double> inverse()``<br>逆行列を求める関数