# 使用するC++コンパイラを指定
CXX = icpc

# コンパイラに渡すフラグ（オプション）
CXXFLAGS = -fopenmp -std=c++17 -Ieigen3 -mkl=sequential -O3

# 生成したい実行可能ファイルの名前
TARGET = occammt

# コンパイル対象のソースファイル
SRC = occammt.cpp

# デフォルトのターゲット
all: $(TARGET)

# 実行可能ファイルを生成するためのルール
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# 'make clean' で生成されたファイルを削除するためのルール
clean:
	rm -f $(TARGET)
