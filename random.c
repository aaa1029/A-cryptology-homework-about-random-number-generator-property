#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// 配置参数
#define TOTAL 5000000    // 总随机数样本量
#define N       4095      // 生成 0 ~ N-1
#define MAX_GAP 20    // 重复间隔最大统计值
#define SEED 42      // 固定随机种子，确保结果可复现

// 全局数据存储
int data[TOTAL];
int freq[N];

// 生成随机数序列：模拟 rand() % N
void generate_data() {
    for (int i = 0; i < TOTAL; i++) {
        data[i] = rand() % N;
        // printf("生成随机数：%d\n", data[i]);
    }
}

// 1. 频数统计
void test_frequency() {
    printf("\n===== 1.1 频数分布统计 =====\n");
    for (int i = 0; i < N; i++) freq[i] = 0;
    for (int i = 0; i < TOTAL; i++) freq[data[i]]++;
    for (int i = 0; i < N-4000; i++) {
        printf("%d → %d 次\n", i, freq[i]);
    }
}

// 2. 概率统计（新增：偏移百分比）
void test_probability() {
    printf("\n===== 1.2 概率分布统计 =====\n");
    double ideal = 1.0 / N;
    printf("理想概率：%.6f\n", ideal);
    for (int i = 0; i < N; i++) {
        double real_p = (double)freq[i] / TOTAL;
        // 计算偏移百分比：(实际概率 - 理想概率) / 理想概率 * 100%
        double offset_percent = (real_p - ideal) / ideal * 100.0;
        
        if (fabs(offset_percent) > 11.0) {  // 只显示偏移超过11%的数值
            printf("%d: %.6f | 偏移: %+.2f%%\n", 
                   i, real_p, offset_percent);
        }
    }
}

// 3. 重复间隔分布
void test_repeat_gap() {
    printf("\n===== 2. 重复间隔分布 =====\n");
    int last[N];
    int gap_count[MAX_GAP] = {0};
    for (int i = 0; i < N; i++) last[i] = -1;

    for (int i = 0; i < TOTAL; i++) {
        int num = data[i];
        if (last[num] != -1) {
            int gap = i - last[num];
            if (gap < MAX_GAP) gap_count[gap]++;
        }
        last[num] = i;
    }

    printf("间隔 1~10 出现次数：\n");
    for (int i = 1; i <= 10; i++) {
        printf("间隔 %d → %d\n", i, gap_count[i]);
    }
}

// 4. 卡方检验
void test_chi_square() {
    printf("\n===== 3. 卡方均匀性检验 =====\n");
    double expect = (double)TOTAL / N;
    double chi = 0;
    for (int i = 0; i < N; i++) {
        chi += pow(freq[i] - expect, 2) / expect;
    }
    printf("期望频数：%.2f\n", expect);
    printf("卡方值：%.2f（越小越均匀）\n", chi);
}

// 5. 游程检验
void test_runs() {
    printf("\n===== 4. 游程检验 =====\n");
    int runs = 1;
    float mean = (N - 1) / 2.0;
    printf("均值：%.2f\n", mean);
    for (int i = 1; i < TOTAL; i++) {
        if ((data[i] > mean)&&(data[i-1] < mean)||(data[i] < mean)&&(data[i-1] > mean)) {
            runs++;
        }
    }
    double expect_runs = TOTAL / 2.0;
    printf("实际游程：%d\n", runs);
    printf("理论游程：%.2f（越接近越随机）\n", expect_runs);
}

// 6. 均值方差检验
void test_mean_variance() {
    printf("\n===== 5. 均值方差检验 =====\n");
    double sum = 0, sum2 = 0;
    for (int i = 0; i < TOTAL; i++) {
        sum += data[i];
        sum2 += data[i] * data[i];
    }
    double mean = sum / TOTAL;
    double var = (sum2 / TOTAL) - (mean * mean);

    double ideal_mean = (N - 1) / 2.0;
    double ideal_var = (N * N - 1) / 12.0;

    printf("实际均值：%.4f | 理论：%.4f\n", mean, ideal_mean);
    printf("实际方差：%.4f | 理论：%.4f\n", var, ideal_var);
}

// 7. 自相关系数检验
void test_autocorr() {
    printf("\n===== 6. 1阶自相关检验 =====\n");
    double mean = 0;
    for (int i = 0; i < TOTAL; i++) mean += data[i];
    mean /= TOTAL;

    double numer = 0, denom = 0;
    for (int i = 0; i < TOTAL - 1; i++) {
        numer += (data[i] - mean) * (data[i+1] - mean);
        denom += (data[i] - mean) * (data[i] - mean);
    }
    double ac = numer / denom;
    printf("自相关系数：%.4f（越接近0越随机）\n", ac);
}

int main() {
    srand(SEED);
    generate_data();

    test_frequency();
    test_probability();
    test_repeat_gap();
    test_chi_square();
    test_runs();
    test_mean_variance();
    test_autocorr();

    printf("\n✅ 6 项测评全部完成！\n");
    return 0;
}