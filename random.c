#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//config
#define TOTAL   5000000 // 总随机数样本量
#define N       10    // 生成 0 ~ N-1
#define MAX_GAP 20      // 重复间隔最大统计值
#define SEED    42     



//全局数据
int data[TOTAL];
int freq[N];
// double normal_data[TOTAL]; // 存储生成的正态分布随机数

void generate_uniform_basic() {
    for (int i = 0; i < TOTAL; i++) {
        data[i] = rand() % N;
        // printf("生成随机数：%d\n", data[i]);
    }
}

void generate_uniform_lcg() {
    // 使用glibc经典LCG参数：周期2^31
    const unsigned int a = 1103515245;
    const unsigned int c = 12345;
    const unsigned int m = 1U << 31;  // 2^31
    unsigned int seed = SEED;

    for (int i = 0; i < TOTAL; i++) {
        seed = (a * seed + c) % m;
        data[i] = seed % N;  // 同样可改用改进版浮点映射
    }
}

// 生成 [0, N-1] 整数的正态分布（中心在 (N-1)/2，拒绝超出范围的样本）
void generate_normal_int() {
    static int has_spare = 0;   // 标记是否有备用正态数
    static double spare = 0.0;  // 保存上次生成的备用正态数

    const double mean = (N - 1) / 2.0;       // 正态分布中心（范围中点）
    const double stddev = (N - 1) / 6.0;     // 标准差：±3σ 覆盖 [0, N-1]

    for (int i = 0; i < TOTAL; i++) {
        int candidate;  // 候选整数
        do {
            double z;   // 标准正态分布数
            if (has_spare) {
                // 直接使用备用正态数
                z = spare;
                has_spare = 0;
            } else {
                double u1, u2, s;
                do {
                    // 生成(0,1)均匀数（避免u1=0导致log(0)错误）
                    u1 = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
                    u2 = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
                    s = u1 * u1 + u2 * u2;
                } while (s >= 1.0);  // 拒绝单位圆外的点

                // Box-Muller变换
                double factor = sqrt(-2.0 * log(s) / s);
                z = u1 * factor;
                spare = u2 * factor;  // 保存另一个正态数备用
                has_spare = 1;
            }

            // 转换为 [0, N-1] 附近的连续值，四舍五入为整数
            double normal_val = mean + stddev * z;
            candidate = (int)round(normal_val);

        } while (candidate < 0 || candidate >= N);  // 拒绝超出范围的样本

        data[i] = candidate;
        // printf("生成正态分布整数：%d\n", data[i]);
    }
}

void test_frequency() {
    // printf("\n1.1 频数分布统计\n");
    for (int i = 0; i < N; i++) freq[i] = 0;
    for (int i = 0; i < TOTAL; i++) freq[data[i]]++;
    for (int i = 0; i < N; i++) {
        // printf("%d : %d 次\n", i, freq[i]);
    }
}

void test_probability() {
    printf("\n1 概率分布统计\n");
    double ideal = 1.0 / N;
    printf("理想概率：%.6f\n", ideal);
    for (int i = 0; i < N; i++) {
        double real_p = (double)freq[i] / TOTAL;
        double offset_percent = ((real_p - ideal) / ideal) * 100.0;
        if (fabs(offset_percent) > 0.0) {  // 只显示偏移超过11%的数值
            printf("%d: %.6f | 偏移: %+.2f%%\n", 
                   i, real_p, offset_percent);
        }
    }
}

void test_repeat_gap() {
    printf("\n2. 重复间隔分布\n");
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
        printf("间隔 %d : %d\n", i, gap_count[i]);
    }
}

void test_chi_square() {
    printf("\n3. 卡方均匀性检验\n");
    double expect = (double)TOTAL / N;
    double chi = 0;
    for (int i = 0; i < N; i++) {
        chi += pow(freq[i] - expect, 2) / expect;
    }
    printf("期望频数：%.2f\n", expect);
    printf("卡方值：%.2f\n", chi);
}

double normal_cdf(double z) {
    return 0.5 * (1.0 + erf(z / sqrt(2.0)));
}


void test_chi_square_normal() {
    printf("\n=== 正态分布卡方检验（修正截断归一化） ===\n");

    const double mean = (N - 1) / 2.0;   // 与生成时一致的均值
    const double stddev = (N - 1) / 6.0; // 与生成时一致的标准差

    // 步骤1：计算截断区间的总概率（归一化用）
    double trunc_low = -0.5;               // 截断左端点（连续性修正）
    double trunc_high = N - 0.5;           // 截断右端点（连续性修正）
    double z_trunc_low = (trunc_low - mean) / stddev;
    double z_trunc_high = (trunc_high - mean) / stddev;
    double p_total = normal_cdf(z_trunc_high) - normal_cdf(z_trunc_low);  // 截断区间总概率

    // 步骤2：计算每个原始组的归一化理论频数
    double *expect = (double *)malloc(N * sizeof(double));
    int *merged = (int *)calloc(N, sizeof(int));
    for (int i = 0; i < N; i++) {
        double z_low = (i - 0.5 - mean) / stddev;   // 区间左端点z-score
        double z_high = (i + 0.5 - mean) / stddev;  // 区间右端点z-score
        double p_original = normal_cdf(z_high) - normal_cdf(z_low);  // 原始概率
        double p_normalized = p_original / p_total;  // 归一化概率（关键修正！）
        expect[i] = p_normalized * TOTAL;  // 归一化后的理论频数
    }

    // 步骤3：合并理论频数 <5 的相邻组（逻辑同前）
    double chi = 0.0;
    int num_groups = 0;
    int k = 0;

    while (k < N) {
        if (merged[k]) { k++; continue; }

        double current_expect = expect[k];
        int current_obs = freq[k];
        int end = k;

        while (current_expect < 5.0 && end + 1 < N) {
            end++;
            merged[end] = 1;
            current_expect += expect[end];
            current_obs += freq[end];
        }

        chi += pow(current_obs - current_expect, 2) / current_expect;
        num_groups++;
        k = end + 1;
    }

    // 步骤4：输出结果
    int df = num_groups - 2 - 1;  // 自由度计算不变
    printf("截断区间：[%.1f, %.1f]，截断总概率：%.4f\n", trunc_low, trunc_high, p_total);
    printf("正态分布参数：均值=%.2f，标准差=%.2f\n", mean, stddev);
    printf("合并后组数：%d，自由度：%d\n", num_groups, df);
    printf("卡方值：%.2f\n", chi);

    if (df > 0) {
        printf("（注：α=0.05, df=%d 时临界值为 %.2f，若卡方值 < 临界值则通过检验）\n", 
               df, 14.07);  // 可根据df替换为对应临界值
    } else {
        printf("警告：自由度≤0，检验无效（建议增加TOTAL）\n");
    }

    free(expect);
    free(merged);
}

void test_runs() {
    printf("\n4. 游程检验 \n");
    int runs = 1;
    float mean = (N - 1) / 2.0;
    printf("均值：%.2f\n", mean);
    for (int i = 1; i < TOTAL; i++) {
        if ((data[i] > mean)&&(data[i-1] < mean)||(data[i] < mean)&&(data[i-1] > mean)) {
            runs++;
        }
    }
    double expect_runs = TOTAL / 2.0;
    printf("理论游程：%.2f\n", expect_runs);
    printf("实际游程：%d\n", runs);
}

void test_mean_and_variance() {
    printf("\n 5. 均值方差检验 \n");
    double sum = 0, sum2 = 0;
    for (int i = 0; i < TOTAL; i++) {
        sum += data[i];
        sum2 += data[i] * data[i];
    }
    double mean = sum / TOTAL;
    double var = (sum2 / TOTAL) - (mean * mean);

    double ideal_mean = (N - 1) / 2.0;
    double ideal_var_uniform = (N * N - 1) / 12.0;
    double ideal_var_normal = pow((N - 1) / 6.0, 2); // 正态分布的方差

    printf("实际均值：%.4f | 理论：%.4f\n", mean, ideal_mean);
    printf("实际方差：%.4f | 理论（均匀分布）：%.4f\n", var, ideal_var_uniform);
    printf("实际方差：%.4f | 理论（正态分布）：%.4f\n", var, ideal_var_normal);
}

void test_autocorr() {
    printf("\n 6. 1阶自相关检验\n");
    double mean = 0;
    for (int i = 0; i < TOTAL; i++) mean += data[i];
    mean /= TOTAL;

    double numer = 0, denom = 0;
    for (int i = 0; i < TOTAL - 1; i++) {
        numer += (data[i] - mean) * (data[i+1] - mean);
    }
    for (int i = 0; i < TOTAL; i++) {
        denom += (data[i] - mean) * (data[i] - mean);
    }
    double ac = numer / denom;
    printf("自相关系数：%.4f\n", ac);
}

void test_all() {
    clock_t start = clock();
    test_frequency();
    test_probability();
    test_repeat_gap();
    test_chi_square();
    test_chi_square_normal();
    test_runs();
    test_mean_and_variance();
    test_autocorr();
    clock_t end = clock();
    printf("\n所有测试耗时: %f 秒\n", ((double)(end - start)) / CLOCKS_PER_SEC);
    return;
}

void generate_data() {
    int flag;
    printf("选择生成方法：1. rand() %% N  2. 线性同余法  3. 正态分布整数\n");
    scanf("%d", &flag);
    clock_t start = clock();
    switch (flag) {
        case 1:
            generate_uniform_basic();
            break;
        case 2:
            generate_uniform_lcg();
            break;
        case 3:
            generate_normal_int();
            break;
        default:
            printf("无效选择，默认使用 rand() %% N\n");
            generate_uniform_basic();
    }
    clock_t end = clock();
    printf("生成数据耗时: %f 秒\n", ((double)(end - start)) / CLOCKS_PER_SEC);
    return;
}

int main(){
    srand(SEED);
    generate_data();
    test_all();
    return 0;
}