# libLinalg

# C++用線形代数ライブラリ

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
