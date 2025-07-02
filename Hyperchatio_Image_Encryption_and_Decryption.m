clc; clear; close all;

%% Load Image
img = imread('C:/Users/soham/Downloads/lena.png');

% Check if the image is RGB and convert to grayscale
if size(img, 3) == 3
    img = rgb2gray(img);
end

img = im2double(img);     
[m, n] = size(img);
num_elements = m * n;

%% Generate Initial Hyperchaotic Sequence
x0 = 0.1234; y0 = 0.5678; z0 = 0.9101; w0 = 0.1121;  
chaotic_seq = hyperchaotic_generator(x0, y0, z0, w0, num_elements);
chaotic_seq = reshape(chaotic_seq, m, n);

%% Generate Quantum-Inspired Dynamic S-Box
SBox = generate_dynamic_sbox(chaotic_seq);

%% Generate Multi-Stage Adaptive Permutation
[~, perm_indices] = sort(chaotic_seq(:));
perm_img = reshape(img(perm_indices), m, n);

%% Bitwise Scrambling
scrambled_img = bitxor(uint8(perm_img * 255), uint8(chaotic_seq * 255));
scrambled_img = double(scrambled_img) / 255;

%% S-Box Substitution (scaled to [0,1])
for i = 1:m
    for j = 1:n
        scrambled_img(i, j) = SBox(mod(floor(scrambled_img(i, j) * 256), 256) + 1) / 255;
    end
end

%% XOR Diffusion + Multi-Round Chaotic Encryption
R = 3;
diffused_img = scrambled_img;
seeds = zeros(R, 4);

for r = 1:R
    x_r = rand; y_r = rand; z_r = rand; w_r = rand;
    seeds(r, :) = [x_r, y_r, z_r, w_r];
    
    chaotic_seq = hyperchaotic_generator(x_r, y_r, z_r, w_r, num_elements);
    chaotic_seq = reshape(chaotic_seq, m, n);

    diffused_img = bitxor(uint8(diffused_img * 255), uint8(chaotic_seq * 255));
    diffused_img = double(diffused_img) / 255;
end

%% Display Encrypted Image
figure;
imshow(img); title('Original Image', 'FontWeight', 'bold');

figure;
imshow(diffused_img); title('Encrypted Image', 'FontWeight', 'bold');

%% ===================== DECRYPTION =====================
enc_img = diffused_img;

%% Reverse XOR Diffusion
for r = R:-1:1
    x_r = seeds(r, 1); y_r = seeds(r, 2);
    z_r = seeds(r, 3); w_r = seeds(r, 4);

    chaotic_seq = hyperchaotic_generator(x_r, y_r, z_r, w_r, num_elements);
    chaotic_seq = reshape(chaotic_seq, m, n);

    enc_img = bitxor(uint8(enc_img * 255), uint8(chaotic_seq * 255));
    enc_img = double(enc_img) / 255;
end

%% Inverse S-Box Substitution
inv_SBox = zeros(1, 256);
for i = 1:256
    inv_SBox(SBox(i) + 1) = i - 1;
end

for i = 1:m
    for j = 1:n
        enc_img(i, j) = inv_SBox(mod(floor(enc_img(i, j) * 256), 256) + 1) / 255;
    end
end

%% Reverse Bitwise Scrambling
chaotic_seq = hyperchaotic_generator(x0, y0, z0, w0, num_elements);
chaotic_seq = reshape(chaotic_seq, m, n);
descrambled_img = bitxor(uint8(enc_img * 255), uint8(chaotic_seq * 255));
descrambled_img = double(descrambled_img) / 255;

%% Reverse Permutation
inv_perm_indices = zeros(size(perm_indices));
inv_perm_indices(perm_indices) = 1:num_elements;
original_img = descrambled_img(inv_perm_indices);
original_img = reshape(original_img, m, n);

%% Display Decrypted Image
figure;
imshow(diffused_img); title('Encrypted Image', 'FontWeight', 'bold');

figure;
imshow(original_img); title('Decrypted Image', 'FontWeight', 'bold');

%% Compute Correlation Coefficients
img_flat = double(img(:));
enc_flat = double(diffused_img(:));
dec_flat = double(original_img(:));

corr_orig_enc = corrcoef(img_flat, enc_flat);
corr_orig_dec = corrcoef(img_flat, dec_flat);
corr_enc_dec  = corrcoef(enc_flat, dec_flat);

fprintf('\nCorrelation (Original vs Encrypted): %.4f\n', corr_orig_enc(1, 2));
fprintf('Correlation (Original vs Decrypted): %.4f\n', corr_orig_dec(1, 2));
fprintf('Correlation (Encrypted vs Decrypted): %.4f\n', corr_enc_dec(1, 2));

%% Plot Histograms
figure;

bins = 0:255;
img_uint8 = uint8(img * 255);
enc_uint8 = uint8(diffused_img * 255);
dec_uint8 = uint8(original_img * 255);

% Original image histogram
subplot(1, 3, 1);
[counts1, ~] = hist(double(img_uint8(:)), bins);
bar(bins, counts1, 'FaceColor', 'b', 'EdgeColor', 'none');
title('Original Image Histogram', 'FontWeight', 'bold');
xlabel('Pixel Intensity'); ylabel('Frequency'); xlim([0 255]);

% Encrypted image histogram
subplot(1, 3, 2);
[counts2, ~] = hist(double(enc_uint8(:)), bins);
bar(bins, counts2, 'FaceColor', 'r', 'EdgeColor', 'none');
title('Encrypted Image Histogram', 'FontWeight', 'bold');
xlabel('Pixel Intensity'); ylabel('Frequency'); xlim([0 255]);

% Decrypted image histogram
subplot(1, 3, 3);
[counts3, ~] = hist(double(dec_uint8(:)), bins);
bar(bins, counts3, 'FaceColor', 'g', 'EdgeColor', 'none');
title('Decrypted Image Histogram', 'FontWeight', 'bold');
xlabel('Pixel Intensity'); ylabel('Frequency'); xlim([0 255]);

%% Histogram Statistics
fprintf('\nHistogram Statistics:\n');
fprintf('Original image - Mean: %.4f, Std Dev: %.4f\n', mean(double(img_uint8(:))), std(double(img_uint8(:))));
fprintf('Encrypted image - Mean: %.4f, Std Dev: %.4f\n', mean(double(enc_uint8(:))), std(double(enc_uint8(:))));
fprintf('Decrypted image - Mean: %.4f, Std Dev: %.4f\n', mean(double(dec_uint8(:))), std(double(dec_uint8(:))));

%% Chi-Square Test & Histogram Intersection
expected = sum(counts2)/length(bins);
chi_square = sum((counts2 - expected).^2 / expected);
fprintf('\nChi-square test for encrypted image uniformity: %.4f\n', chi_square);

hist_intersection = sum(min(counts1/sum(counts1), counts3/sum(counts3)));
fprintf('Histogram intersection (Original vs Decrypted): %.4f\n', hist_intersection);

%% ============ Helper Functions ============

function seq = hyperchaotic_generator(x0, y0, z0, w0, N)
    seq = zeros(1, N);
    a = 35; b = 3; c = 28; d = 5;
    dt = 0.01;
    for i = 1:N
        x0 = x0 + dt * (a * (y0 - x0));
        y0 = y0 + dt * (b * x0 - y0 - x0 * z0);
        z0 = z0 + dt * (c * x0 * y0 - d * z0);
        w0 = w0 + dt * (x0 * y0 - w0);
        seq(i) = mod(x0 + y0 + z0 + w0, 1);
    end
end

function SBox = generate_dynamic_sbox(chaotic_seq)
    chaotic_seq = mod(floor(chaotic_seq(:) * 256), 256);
    [~, indices] = sort(chaotic_seq(1:256));
    unique_vals = unique(indices, 'stable');

    if length(unique_vals) < 256
        missing = setdiff(0:255, unique_vals);
        unique_vals = [unique_vals(:); missing(:)];
    elseif length(unique_vals) > 256
        unique_vals = unique_vals(1:256);
    end

    SBox = zeros(1, 256);
    SBox(unique_vals + 1) = 0:255;

    if length(unique(SBox)) < 256
        error('Generated SBox is not bijective.');
    end
end
