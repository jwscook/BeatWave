NT = 2^9


using Plots, FFTW, PyCall

sdf = pyimport("sdf")
F = zeros(8192, NT);

for i in 0:size(F, 2)-1
  h = sdf.read(lpad(i, 5, "0") * ".sdf")
  F[:, i+1] .= h.Electric_Field_Ex_averaged.data
end

filt = sin.(((1:size(F,2)) .- 0.5) ./ size(F,2) .* pi)';
heatmap(F[1:end÷8, :])

# heatmap(log10.(abs.((fft(F .* filt))[1:32, 1:end÷2]))')



