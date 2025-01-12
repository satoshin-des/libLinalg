# libLinalg

# C++用線形代数ライブラリ

## 全体的に準備中

## 関数について

### Vector classの関数
- ``print``：ベクトルを表示する関数
- ``setLength(const int length)``：ベクトルの長さを``length``に設定する関数
- ``int length()``：ベクトルの長さを返す関数
- ``Vector<typename> cast<typename>()``：ベクトルをcastする関数
- ``at(const int index)``：ベクトルの``index``番目の要素を取ってくる関数
- ``bool isZero()``：ベクトルが零ベクトルであるかを判定する関数
- ``innerProduct(Vector<T> vector)``：``vector``との標準内積を求める関数
- ``norm(const unsigned int p = 2)``：$l_p$-ノルムを求める関数