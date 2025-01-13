# C++用線形代数ライブラリ「libLinalg」

## 全体的に準備中

## 関数について

### Vector classの関数

- ``+``，``-``
それぞれ，ベクトル同士の和と差
- ``*``
ベクトルとスカラーの場合はスカラー倍，ベクトルとベクトルの場合は内積
- ``/``
ベクトルの各要素を実数で割る
- ``%``
ベクトルの各要素を整数で割った余り
- ``void print()``
ベクトルを表示する関数
- ``setLength(const int length)``
ベクトルの長さを``length``に設定する関数
- ``int length()``：ベクトルの長さを返す関数
- ``Vector<U> cast<U>()``
ベクトルをcastする関数
- ``substitute(const int index, const T number_to_substitute)``
ベクトルの``index``番目に``number_to_substitute``を代入する関数
- ``substitute(std::vector<T> vector)``
ベクトルに``vector``を代入する関数
- ``int at(const int index)``
ベクトルの``index``番目の要素を取ってくる関数
- ``bool isZero()``
ベクトルが零ベクトルであるかを判定する関数
- ``void setSeed(const unsigned int seed)``
乱数のシードを設定する関数
- ``void random()``
ベクトルをランダムに生成する関数
- ``void random(const double lower_bound, const double upper_bound)``
ベクトルをランダムに生成する関数．但し，各要素は``lower_bound``以上``upper_bound``以下である
- ``T innerProduct(Vector<T> vector)``
``vector``との標準内積を求める関数
- ``norm(const unsigned int p = 2)``
$l_p$-ノルムを求める関数

### Matrix classの関数
- ``Matrix<T> power(const int n)``
行列の``n``乗を計算する関数
- ``void print()``
行列を表示する関数
- ``T at(const int row_index, const int column_index)``
行列の(``row_index``, ``column_index``)-成分を返す関数
- ``void setSeed(const unsigned int seed)``
乱数のシードを設定する関数
- ``void random()``
行列をランダムに生成する関数
- ``void random(const double lower_bound, const double upper_bound)``
行列をランダムに生成する関数．但し，各要素は``lower_bound``以上``upper_bound``以下である
- ``void setDimensions(const int num_rows, num int num_columns)``
行列のサイズを``num_rows``×``num_columns``に設定する関数
- ``void setIdentity(const int size)``
行列を``size``次単位行列に設定する関数
- ``bool isRegular()``
行列が正則かどうか判定する関数
- ``bool isSquare()``
行列が正方かどうか判定する関数
- ``void substitute(const std::vector<std::vector<T>> matrix)``
行列に``matrix``を代入する関数
- ``void substituteRow(const int index, std::vector<T> vector)``
- ``void substituteRow(const int index, Vector<T> vector)``
行列の``index``行目にベクトル``vector``を代入する関数
- ``void substituteColumn(const int index, std::vector<T> vector)``
- ``void substituteColumn(const int index, Vector<T> vector)``
行列の``index``列目にベクトル``vector``を代入する関数
- ``void substitute(const int row_index, const int column_index, const T number)``
行列の(``row_index``，``column_index``)-成分に``number``を代入する関数
- ``Matrix<U> cast<U>()``
行列をcastする関数
- ``Vector<T> row(const int index)``
行列の第``index``行を取り出す関数
- ``int rows()`
行列の行数を返す関数
- ``int columns()``
行列の列数を返す関数
- ``T det()``
- ``T determinant()``
行列の行列式を求める関数
- ``Matrix<T> transpose()``
- ``Matrix<T> trans()``
転置行列を求める関数
- ``T trace()``
行列のトレースを求める関数
- ``T coFactor(const int row_index, const int column_index)``
行列の(``row_index``, ``column_index``)余因子を求める関数
- ``Matrix<T> adjugate()``
行列の余因子行列を求める関数
- ``Matrix<double> inverse()``
逆行列を求める関数