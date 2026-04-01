# Analisis Numerik Nilai Eigen Matriks

Proyek ini berisi implementasi Octave/MATLAB untuk menghitung nilai eigen dan vektor eigen matriks persegi dengan pendekatan:

- bisection pada fungsi karakteristik
- verifikasi turunan numerik untuk mendeteksi akar kembar
- eliminasi Gauss pivoting untuk konstruksi vektor eigen

## Latar Belakang Singkat

Masalah nilai eigen ditulis sebagai:

$$
Av = \lambda v \quad \Rightarrow \quad (A - \lambda I)v = 0
$$

Nilai eigen diperoleh dari akar-akar:

$$
f(\lambda) = \det(A - \lambda I) = 0
$$

## Ide Metode

1. Tentukan domain pencarian nilai eigen dari diagonal dan simpangan baris:

$$
R_i = \sum_{j \ne i} |a_{ij}|
$$

$$
\lambda_{\min} = \min_i(a_{ii} - R_i) - 1, \qquad
\lambda_{\max} = \max_i(a_{ii} + R_i) + 1
$$

2. Diskritisasi domain menjadi segmen kecil.
3. Deteksi akar sederhana dari syarat:

$$
f(a)f(b) \le 0
$$

4. Untuk akar kembar, cari kandidat dari turunan numerik:

$$
f'(\lambda) \approx \frac{f(\lambda + h) - f(\lambda)}{h}
$$

5. Verifikasi kandidat akar kembar dengan:

$$
|f(\lambda_{ext})| < 10^{-4}
$$

## Struktur File

- `eigen.m` : fungsi utama `MyEigen` + helper internal
- `MyBisection.m` : metode bisection
- `MyDet.m` : determinan via LU
- `MyLU.m` : dekomposisi LU dengan pivoting
- `MyLUX.m` : solver sistem linear dari LU
- `GaussPivot.m` : eliminasi Gauss pivoting
- `main.tex` : laporan LaTeX
- `main.pdf` : hasil kompilasi laporan

## Cara Menjalankan (Octave/MATLAB)

```matlab
clc; clear;

A = [3 0 0;
     1 4 -2;
    -1 -1 5];

[V, D] = MyEigen(A, 1e-8, 200);

disp('Matriks eigenvector V:');
disp(V);

disp('Matriks eigenvalue D:');
disp(D);
```

## Contoh Hasil Uji

Untuk matriks contoh di atas, eigenvalue teoritis adalah:

$$
\lambda = 3,\;3,\;6
$$

Sehingga output diagonal matriks $D$ diharapkan mendekati:

$$
D = \mathrm{diag}(3, 3, 6)
$$

## Catatan

- Pendekatan ini fokus pada kestabilan numerik dan deteksi akar kembar pada kasus nilai eigen real.
- Untuk akar kembar, kolom vektor eigen bisa tampak identik pada konstruksi tertentu, dan ini dapat terjadi secara numerik tergantung struktur matriks.
