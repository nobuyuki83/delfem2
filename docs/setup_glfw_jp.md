# GLFWライブラリの設定方法



GLFWは以下の３つの方法で設定できます．

- パッケージ・マネージャを使って `glfw` を設定する (MacかUbuntu向け)
- 予めコンパイルされたバイナリをダウンロードする(MacかWindows向け)
- GLFWのコードを自分でコンパイルする

以下にそれぞれの方法について説明します．

----



##  パッケージ・マネージャを使う方法

MacOSでは，`brew`というパッケージマネージャを使って`glfw`をインストールできます．

```bash
$ brew install glfw
```
Ubuntuでは`apt-get`というパッケージマネージャを使って `glfw` をインストールできます．

```bash
$ sudo apt-get install -y libx11-dev xorg-dev \
                          libglu1-mesa libglu1-mesa-dev \
                          libgl1-mesa-glx libgl1-mesa-dev
$ sudo apt install -y libglfw3 libglfw3-dev
$ sudo apt install -y libglew-dev
```
Windowsにはこのようなパッケージマネージャはなさそうです．

---



## 予めコンパイルされたバイナリをダウンロード

WindowやMac向けにコンパイル済みのバイナリライブラリが提供されています．

 https://www.glfw.org/download.html

圧縮されたファイルを展開して， `GLFW_Lib` というフォルダに名前を変更し， `3rd_party/`フォルダの下において下さい．

ヘッダファイル `glfw3.h` が以下の場所に置かれるはずです

```
ProjectTemplate_CppGlfw/3rd_party/GLFW_Lib/include/GLFW/glfw3.h
```





---


## ソースコードを自分でビルド

以下のコマンドで，ソースコードから `glfw` ライブラリをコンパイルしてフォルダ `3rd_party/GLFW_Lib` の中に置くことができます．

```bash
$ mkdir 3rd_party/GLFW_Lib 
$ git submodule update --init 3rd_party/glfw
$ cd 3rd_party/glfw
$ cmake .
$ cmake --build . --config Release
$ cmake --install . --prefix ../GLFW_Lib 
```

ヘッダファイル `glfw3.h` が以下の場所に置かれるはずです

```
ProjectTemplate_CppGlfw/3rd_party/GLFW_Lib/include/GLFW/glfw3.h
```







 



