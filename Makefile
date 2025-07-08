# 使用するC++コンパイラを指定
CXX = g++

# コンパイラに渡すフラグ（オプション）
CXXFLAGS = -std=c++17 -I./eigen3 -I./matplotlib-cpp $(shell python3.10-config --cflags)

# リンカに渡すフラグ（オプション）
# --ldflagsの出力に-lpython3.10が不足しているため、手動で追加
LDFLAGS = $(shell python3.10-config --ldflags) -lpython3.10 -lpthread

# 生成したい実行可能ファイルの名前
TARGET = occam15d

# コンパイル対象のソースファイル
SRC = occam15d.cpp

# デフォルトのターゲット
all: $(TARGET)

# 実行可能ファイルを生成するためのルール
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# 'make clean' で生成されたファイルを削除するためのルール
clean:
	rm -f $(TARGET)
