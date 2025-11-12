/*
    円盤の積み上げ
    落下速度と壁の速度をパラメータにする
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>

#define PI atan(1.0) * 4
#define ratio 10.0

/*=========物性パラメータ（無次元化）==========*/
#define N 3             // 粒子の数
#define wall 2          // 壁の数
#define r 1.0           // radius
#define e 0.3           // coefficient of restitution
#define mu 0.7          // the friction coefficient
#define kn 1e+2 * ratio // the normal elastic constant
#define kt 4e+1 * ratio // the tangential elastic constant
#define g 1.0           // gravity
#define vfall 0.1       /*粒子を落とす速度*/
#define hight 0.1       // 円盤2を落とす高さ

/*===========壁のパラメータ=============*/
#define A 1.2e-4   /*押し込みの幅 k=1e+2 A=1e-3, k=1e+3 A=1.2e-4, k=1e+4 A=2e-5*/
#define vwall 1e-6 /*壁の開放速度*/

/*============時間のパラメータ=============*/
// #define dt sqrt(1 / kn) / 1000 // 時間刻み幅
#define dt 1e-3     // 時間刻み幅
#define data 200    /*gif画像の枚数*/
#define tend 600000 /*シミュレーション時間*/

void init(int n, double *a, double *b, double *c, double x, double y, double z)
{
    a[n] = x;
    b[n] = y;
    c[n] = z;
};
void vec_init0(int n, double o[n]) // n行ベクトルoを0で初期化
{
    int i = 0;
    for (i = 0; i < n; i++)
    {
        o[i] = 0;
    }
};
void matrix_init0(int m, int n, double o[m][n]) // m×n行列oを0で初期化
{
    int i = 0, j = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            o[i][j] = 0;
        }
    }
};

int main(void)
{
    FILE *fp, *fp1, *fp2, *fp3;
    int count = 1, i, j, k, igraph, jgraph, kgraph;
    double t, tstar; // 時刻
    int it;          // ステップ数
    double tau;      // 接触する時刻

    /*==============粒子座標の定義==============*/
    double x[N], xnew[N], y[N], ynew[N], theta[N], thetanew[N]; // 粒子の位置
    double u[N], unew[N], v[N], vnew[N], omega[N], omeganew[N]; // 粒子の速度

    double fx[N], fy[N];

    /*===========粒子間相互作用の定義============*/
    double l[N][N];                                // ２粒子間の距離
    double lx[N][N], ly[N][N], du[N][N], dv[N][N]; // 相対距離と相対速度
    double delta_nx[N][N], delta_ny[N][N];         // 法線方向相対変位ベクトル
    double v_nx[N][N], v_ny[N][N];                 // 法線方向相対速度ベクトル
    double nx[N][N], ny[N][N];                     // 法線ベクトル
    double fnx[N][N], fny[N][N];                   // 法線方向合力
    double delta_tx[N][N], delta_ty[N][N];
    double delta_txnew[N][N], delta_tynew[N][N];
    double deltat[N][N];
    double vtx[N][N], vty[N][N];
    double vt[N][N];
    double tx[N][N], ty[N][N];
    double ftx[N][N], fty[N][N];
    double T[N];
    double ft[N][N], fn[N][N];

    /*===============壁との相互作用==============*/
    double delta_nx_wall[N][wall], delta_ny_wall[N][wall];
    double vnx_wall[N][wall], vny_wall[N][wall];
    double fnx_wall[N][wall], fny_wall[N][wall];
    double l_wall[N][wall];
    double delta_tx_wall[N][wall], delta_ty_wall[N][wall];
    double delta_txnew_wall[N][wall], delta_tynew_wall[N][wall];
    double deltat_wall[N][wall];
    double vtx_wall[N][wall], vty_wall[N][wall];
    double vt_wall[N][wall];
    double tx_wall[N][wall], ty_wall[N][wall];
    double ftx_wall[N][wall], fty_wall[N][wall];
    double T_wall[N][wall];
    double ft_wall[N][wall], fn_wall[N][wall];

    /*===============壁の定義===================*/
    double nx_wall[wall], ny_wall[wall], d[wall]; // 法線ベクトルと切片
    double u_wall[wall], v_wall[wall];            // 壁の速度

    /*=============物性パラメータ=================*/
    double m = 1.0;                                                         // 質量
    double I = 0.5 * m * r * r;                                             // 慣性モーメント
    double etan = -2 * log(e) * sqrt(m * kn / (PI * PI + log(e) * log(e))); // 粘性係数(法線方向)
    double etat = -2 * log(e) * sqrt(m * kt / (PI * PI + log(e) * log(e))); // 粘性係数(接線方向)
    double mu_wall = 0.2;

    /*==========時間刻み，プロットする点の間隔=======*/
    int itend = tend / dt;
    int countgraph = (int)(tend / dt / data); /*countgraphごとにファイルへ出力*/
    /*=============ファイルの指定==================*/
    char filename[100];
    sprintf(filename, "data/k%.1e.txt", kn);
    fp = fopen(filename, "w");

    for (mu_wall = 0.1; mu_wall < 0.25; mu_wall += 0.01)
    {
        int flag = 0;

        /*降伏応力*/
        double fmax = -1;

        /*================初期条件===================*/
        t = 0.0;
        it = 0;

        FILE *file;

        /*======位置と速度の読み込み=======*/
        sprintf(filename, "init/init_p_A%.1e.txt", A);
        file = fopen(filename, "r");
        // ファイルから値を読み込む
        for (i = 0; i < N; i++)
        {
            fscanf(file, "%lf %lf %lf %lf %lf %lf", &x[i], &y[i], &theta[i], &u[i], &v[i], &omega[i]);
        }

        // ファイルを閉じる
        fclose(file);

        /*======delta_tの読み込み=======*/
        sprintf(filename, "init/init_deltat_A%.1e.txt", A);
        file = fopen(filename, "r");

        if (file == NULL)
        {
            printf("ファイルが開けません。\n");
            return 1;
        }

        // ファイルから値を読み込む
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                fscanf(file, "%lf", &delta_tx[i][j]);
            }
        }

        // ファイルから値を読み込む
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                fscanf(file, "%lf", &delta_ty[i][j]);
            }
        }
        for (i = 0; i < N; i++)
        {

            fscanf(file, "%lf", &delta_tx_wall[i][0]);
        }
        for (i = 0; i < N; i++)
        {

            fscanf(file, "%lf", &delta_ty_wall[i][0]);
        }

        // ファイルを閉じる
        fclose(file);

        matrix_init0(N, N, delta_txnew);
        matrix_init0(N, N, delta_tynew);

        /*壁の初期条件*/
        init(0, nx_wall, ny_wall, d, 0.0, 1.0, 0.0);
        u_wall[0] = 0.0;
        v_wall[0] = 0.0;
        init(1, nx_wall, ny_wall, d, 0.0, -1.0, 2 + sqrt(3));
        u_wall[1] = 0.0;
        v_wall[1] = 0.0;
        matrix_init0(N, wall, delta_txnew_wall);
        matrix_init0(N, wall, delta_tynew_wall);

        /*==================時間発展==================*/
        while (flag == 0)
        {
            if (t < 1)
            {
                d[1] = 2 + sqrt(3);
                v_wall[1] = 0.0;
            }
            else
            {
                d[1] = 2 + sqrt(3) - vwall * (t - 1);
                v_wall[1] = -vwall;
            }

            for (i = 0; i < N; i++)
            {
                /*============粒子iに働く力の相互作用計算=============*/
                for (j = 0; j < N; j++) // 粒子jとの相互作用
                {
                    if (j == i)
                    {
                        continue;
                    }
                    // printf("%d %d\n", i, j);

                    /*==================距離と三角関数===================*/
                    lx[i][j] = x[i] - x[j];
                    ly[i][j] = y[i] - y[j];
                    du[i][j] = u[i] - u[j];
                    dv[i][j] = v[i] - v[j];
                    l[i][j] = sqrt(pow(lx[i][j], 2) + pow(ly[i][j], 2));

                    /*=====================法線ベクトル==================*/
                    nx[i][j] = lx[i][j] / l[i][j];
                    ny[i][j] = ly[i][j] / l[i][j];

                    /*===============法線方向変位ベクトル=================*/
                    delta_nx[i][j] = (l[i][j] - 2 * r) * nx[i][j];
                    delta_ny[i][j] = (l[i][j] - 2 * r) * ny[i][j];

                    /*================法線方向相対速度ベクトル=============*/
                    v_nx[i][j] = (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * nx[i][j];
                    v_ny[i][j] = (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * ny[i][j];

                    /*================接線方向相対速度ベクトル=============*/
                    vtx[i][j] = du[i][j] - (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * nx[i][j] + r * (omega[i] + omega[j]) * ny[i][j];
                    vty[i][j] = dv[i][j] - (du[i][j] * nx[i][j] + dv[i][j] * ny[i][j]) * ny[i][j] - r * (omega[i] + omega[j]) * nx[i][j];

                    /*==============相対変位ベクトル大きさ=================*/
                    deltat[i][j] = sqrt(pow(delta_tx[i][j], 2) + pow(delta_ty[i][j], 2));
                    vt[i][j] = sqrt(pow(vtx[i][j], 2) + pow(vty[i][j], 2));

                    /*================接線ベクトル=======================*/
                    if (vt[i][j] == 0 && deltat[i][j] > 0)
                    {
                        tx[i][j] = delta_tx[i][j] / deltat[i][j];
                        ty[i][j] = delta_ty[i][j] / deltat[i][j];
                    }
                    else if (vt[i][j] > 0)
                    {
                        tx[i][j] = vtx[i][j] / vt[i][j];
                        ty[i][j] = vty[i][j] / vt[i][j];
                    }
                    else
                    {
                        tx[i][j] = 0;
                        ty[i][j] = 0;
                    }

                    /*==================法線方向接触力の更新================*/
                    if (l[i][j] - 2 * r < 0) // 接触
                    {
                        fnx[i][j] = -kn * delta_nx[i][j] - etan * v_nx[i][j];
                        fny[i][j] = -kn * delta_ny[i][j] - etan * v_ny[i][j];
                    }
                    else // 非接触
                    {
                        fnx[i][j] = 0;
                        fny[i][j] = 0;
                    }

                    /*===============接線方向接触力の予測====================*/
                    if (l[i][j] - 2 * r < 0)
                    {
                        ftx[i][j] = -kt * delta_tx[i][j] - etat * vtx[i][j];
                        fty[i][j] = -kt * delta_ty[i][j] - etat * vty[i][j];
                    }
                    else
                    {
                        ftx[i][j] = 0;
                        fty[i][j] = 0;
                    }

                    ft[i][j] = sqrt(pow(ftx[i][j], 2) + pow(fty[i][j], 2));
                    fn[i][j] = sqrt(pow(fnx[i][j], 2) + pow(fny[i][j], 2));

                    /*============接線方向接触力の決定,累積滑り変位の更新========*/
                    if (ft[i][j] >= mu * fn[i][j] && l[i][j] - 2 * r < 0) // slip かつ 接触
                    {
                        ftx[i][j] = -mu * fn[i][j] * tx[i][j];
                        fty[i][j] = -mu * fn[i][j] * ty[i][j];

                        /*============deltatを上限値で打ち切る==============*/
                        delta_txnew[i][j] = mu * kn * (2 * r - l[i][j]) / kt * tx[i][j];
                        delta_tynew[i][j] = mu * kn * (2 * r - l[i][j]) / kt * ty[i][j];
                    }
                    else if (ft[i][j] < mu * fn[i][j] && l[i][j] - 2 * r < 0) // no-slip かつ 接触
                    {
                        delta_txnew[i][j] = delta_tx[i][j] + vtx[i][j] * dt;
                        delta_tynew[i][j] = delta_ty[i][j] + vty[i][j] * dt;
                    }
                    else // 非接触
                    {
                        delta_tx[i][j] = 0;
                        delta_ty[i][j] = 0;
                    }

                    delta_tx[i][j] = delta_txnew[i][j];
                    delta_ty[i][j] = delta_tynew[i][j];
                }

                /*==================壁との相互作用====================*/
                for (k = 0; k < wall; k++)
                {
                    /*壁kとの距離*/
                    l_wall[i][k] = nx_wall[k] * x[i] + ny_wall[k] * y[i] + d[k];

                    /*法線方向変位ベクトル*/
                    delta_nx_wall[i][k] = (l_wall[i][k] - r) * nx_wall[k];
                    delta_ny_wall[i][k] = (l_wall[i][k] - r) * ny_wall[k];

                    /*法線方向相対速度ベクトル*/
                    vnx_wall[i][k] = ((u[i] - u_wall[k]) * nx_wall[k] + (v[i] - v_wall[k]) * ny_wall[k]) * nx_wall[k];
                    vny_wall[i][k] = ((u[i] - u_wall[k]) * nx_wall[k] + (v[i] - v_wall[k]) * ny_wall[k]) * ny_wall[k];

                    /*接線方向相対ベクトル*/
                    vtx_wall[i][k] = (u[i] - u_wall[k]) - vnx_wall[i][k] + r * omega[i] * ny_wall[k];
                    vty_wall[i][k] = (v[i] - v_wall[k]) - vny_wall[i][k] - r * omega[i] * nx_wall[k];

                    /*相対変位ベクトル*/
                    deltat_wall[i][k] = sqrt(pow(delta_tx_wall[i][k], 2) + pow(delta_ty_wall[i][k], 2));
                    vt_wall[i][k] = sqrt(pow(vtx_wall[i][k], 2) + pow(vty_wall[i][k], 2));

                    /*接線ベクトル*/
                    if (vt_wall[i][k] == 0 && deltat_wall[i][k] > 0)
                    {
                        tx_wall[i][k] = delta_tx_wall[i][k] / deltat_wall[i][k];
                        ty_wall[i][k] = delta_ty_wall[i][k] / deltat_wall[i][k];
                    }
                    else if (vt_wall[i][k] > 0)
                    {
                        tx_wall[i][k] = vtx_wall[i][k] / vt_wall[i][k];
                        ty_wall[i][k] = vty_wall[i][k] / vt_wall[i][k];
                    }
                    else
                    {
                        tx_wall[i][k] = 0;
                        ty_wall[i][k] = 0;
                    }

                    /*法線方向接触力の更新*/
                    if (l_wall[i][k] - r < 0)
                    {
                        fnx_wall[i][k] = -kn * delta_nx_wall[i][k] - etan * vnx_wall[i][k];
                        fny_wall[i][k] = -kn * delta_ny_wall[i][k] - etan * vny_wall[i][k];
                    }
                    else
                    {
                        fnx_wall[i][k] = 0;
                        fny_wall[i][k] = 0;
                    }

                    /*接線方向力(摩擦力)の予測*/
                    if (l_wall[i][k] - r < 0)
                    {
                        ftx_wall[i][k] = -kt * delta_tx_wall[i][k] - etat * vtx_wall[i][k];
                        fty_wall[i][k] = -kt * delta_ty_wall[i][k] - etat * vty_wall[i][k];
                    }
                    else
                    {
                        ftx_wall[i][k] = 0;
                        fty_wall[i][k] = 0;
                    }

                    fn_wall[i][k] = sqrt(pow(fnx_wall[i][k], 2) + pow(fny_wall[i][k], 2));
                    ft_wall[i][k] = sqrt(pow(ftx_wall[i][k], 2) + pow(fty_wall[i][k], 2));

                    /*接線方向接触力の決定、累積すべり変位の更新*/
                    if (ft_wall[i][k] >= mu_wall * fn_wall[i][k] && l_wall[i][k] - r < 0) // slip かつ 接触
                    {
                        ftx_wall[i][k] = -mu_wall * fn_wall[i][k] * tx_wall[i][k];
                        fty_wall[i][k] = -mu_wall * fn_wall[i][k] * ty_wall[i][k];

                        // printf("%f %f\n", t, ftx_wall[i][k]);

                        /*============deltatを上限値で打ち切る==============*/
                        delta_txnew_wall[i][k] = mu_wall * kn * (r - l_wall[i][k]) / kt * tx_wall[i][k];
                        delta_tynew_wall[i][k] = mu_wall * kn * (r - l_wall[i][k]) / kt * ty_wall[i][k];
                    }
                    else if (ft_wall[i][k] < mu_wall * fn_wall[i][k] && l_wall[i][k] - r < 0) // no-slip かつ 接触
                    {
                        delta_txnew_wall[i][k] = delta_tx_wall[i][k] + vtx_wall[i][k] * dt;
                        delta_tynew_wall[i][k] = delta_ty_wall[i][k] + vty_wall[i][k] * dt;
                    }
                    else
                    {
                        delta_txnew_wall[i][k] = 0;
                        delta_tynew_wall[i][k] = 0;
                    }

                    delta_tx_wall[i][k] = delta_txnew_wall[i][k];
                    delta_ty_wall[i][k] = delta_tynew_wall[i][k];
                }
            }

            /*===========位置と速度の更新　ファイルへの記録==============*/
            for (i = 0; i < N; i++)
            {
                fx[i] = 0.0;
                fy[i] = 0.0;
                T[i] = 0.0;
            }

            for (i = 0; i < N; i++)
            {

                /*粒子間相互作用を合力に加算*/
                for (j = 0; j < N; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    fx[i] += fnx[i][j] + ftx[i][j];
                    fy[i] += fny[i][j] + fty[i][j];
                    T[i] += -r * (nx[i][j] * fty[i][j] - ny[i][j] * ftx[i][j]);
                }

                /*壁との相互作用を合力に加算*/
                for (k = 0; k < wall; k++)
                {
                    fx[i] += fnx_wall[i][k] + ftx_wall[i][k];
                    fy[i] += fny_wall[i][k] + fty_wall[i][k];
                    T[i] += -r * (nx_wall[k] * fty_wall[i][k] - ny_wall[k] * ftx_wall[i][k]);
                }

                /*重力を合力に加算*/
                fy[i] += -g;

                /*=========================位置と速度の更新======================*/

                /*速度更新(一次のシンプレクティックorリープフロッグ)*/
                unew[i] = u[i] + fx[i] * dt;
                vnew[i] = v[i] + fy[i] * dt;
                omeganew[i] = omega[i] + T[i] * dt / I;

                /*位置更新(一次のシンプレクティックorリープフロッグ)*/
                xnew[i] = x[i] + unew[i] * dt;
                ynew[i] = y[i] + vnew[i] * dt;
                thetanew[i] = theta[i] + omeganew[i] * dt;

                x[i] = xnew[i];
                y[i] = ynew[i];
                u[i] = unew[i];
                v[i] = vnew[i];
                theta[i] = thetanew[i];
                omega[i] = omeganew[i];
            }

            if (fmax > -fny_wall[2][1])
            {
                flag = 1;
            }

            /*応答の最大値*/
            if (fmax < -fny_wall[2][1])
            {
                fmax = -fny_wall[2][1];
            }

            /*===========時間の更新==========*/
            t += dt;
            count++;
        }

        fprintf(fp, "%f %f\n", mu_wall, fmax);
    }

    return 0;
}