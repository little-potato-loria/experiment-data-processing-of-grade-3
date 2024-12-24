#include <stdio.h>
#include <stdlib.h>

// 比较两个字节并返回不同的位数
int count_different_bits(unsigned char byte1, unsigned char byte2) {
    unsigned char diff = byte1 ^ byte2;  // 异或运算得到不同的位
    int count = 0;
    while (diff) {
        count += diff & 1;  // 计算低位的1
        diff >>= 1;  // 右移一位
    }
    return count;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "用法: %s file1 file2\n", argv[0]);
        return 1;
    }

    const char *file1 = argv[1];
    const char *file2 = argv[2];
    FILE *f1 = fopen(file1, "rb");
    FILE *f2 = fopen(file2, "rb");

    if (!f1 || !f2) {
        fprintf(stderr, "无法打开其中一个文件\n");
        if (f1) fclose(f1);
        if (f2) fclose(f2);
        return 1;
    }

    int diff_bit_count = 0;
    int total_bits = 0;
    int ch1, ch2;
    while ((ch1 = fgetc(f1)) != EOF && (ch2 = fgetc(f2)) != EOF) {
        diff_bit_count += count_different_bits(ch1, ch2);
        total_bits += 8;  // 每个字节有8位
    }

    if (fgetc(f1) != EOF || fgetc(f2) != EOF) {
        fprintf(stderr, "文件长度不同\n");
    } else {
        printf("文件 %s 和 %s 有 %d 个不同的位\n", file1, file2, diff_bit_count);
        printf("文件总共有 %d 位\n", total_bits);
    }

    fclose(f1);
    fclose(f2);

    return 0;
}

