#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TOTAL 100  // 总随机数样本量
#define N 300          // 生成 0~N-1

// ===================== 工具1：频数统计 =====================
void test_frequency(int *freq) {
    printf("\n===== 工具1：频数分布统计 =====\n");
    for (int i = 0; i < N; i++) {
        printf("%d → %d 次\n", i, freq[i]);
    }
}

// ===================== 工具2：概率统计 =====================
void test_probability(int *freq) {
    printf("\n===== 工具2：概率分布统计 =====\n");
    for (int i = 0; i < N; i++) {
        double p = (double)freq[i] / TOTAL;
        printf("%d → 概率 %.6f\n", i, p);
    }
}

// ===================== 工具3：重复间隔分布 =====================
void test_repeat_gap() {
    printf("\n===== 工具3：重复间隔分布 =====\n");
    int last[N] = {-1};
    int gap_count[1000] = {0};  // 间隔统计
    int prev = rand() % N;
    last[prev] = 0;

    for (int i = 1; i < TOTAL; i++) {
        int num = rand() % N;
        if (last[num] != -1) {
            int gap = i - last[num];
            if (gap < 1000) gap_count[gap]++;
        }
        last[num] = i;
    }

    printf("间隔 1~10 的出现次数：\n");
    for (int i = 1; i <= 10; i++) {
        printf("间隔 %d → %d 次\n", i, gap_count[i]);
    }
}

// ===================== 工具4：卡方检验 =====================
void test_chi_square(int *freq) {
    printf("\n===== 工具4：卡方均匀性检验 =====\n");
    double expect = (double)TOTAL / N;
    double chi = 0;
    for (int i = 0; i < N; i++) {
        chi += pow(freq[i] - expect, 2) / expect;
    }
    printf("卡方值 = %.2f\n", chi);
    printf("越小越均匀，越大越不均匀\n");
}

// ===================== 工具5：游程检验 =====================
void test_runs() {
    printf("\n===== 工具5：游程检验（独立性） =====\n");
    int prev = rand() % N;
    int runs = 1;
    for (int i = 1; i < TOTAL; i++) {
        int curr = rand() % N;
        if ((curr > prev && prev <= curr) || (curr < prev && prev >= curr)) {
            runs++;
        }
        prev = curr;
    }
    printf("总游程数 = %d\n", runs);
    double expect_runs = (2.0 * TOTAL - 1) / 3.0;
    printf("理论期望游程 = %.2f\n", expect_runs);
    printf("接近则随机好，偏差大则有规律\n");
}

// ===================== 主函数 =====================
int main() {
    srand((unsigned)time(NULL));
    int freq[N] = {0};

    // 生成大量随机数
    for (int i = 0; i < TOTAL; i++) {
        int num = rand() % N;
        printf("生成随机数：%d\n", num);
        freq[num]++;
    }

    // 5大测评工具调用
    test_frequency(freq);
    test_probability(freq);
    test_repeat_gap();
    test_chi_square(freq);
    test_runs();

    return 0;
}