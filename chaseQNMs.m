close
clear all
format longE

addpath('/path/to/advanpix/');

list = [200, 150, 100];
L = length(list);

Marks = {'r.', 'b*', 'kd'};

if length(Marks) < L
  error('Error: Vector of markers is insufficient.');
end % if

qnms = cell(1, L);
legs = cell(1, L);

for idx = 1:L
  n = list(idx);
  nstr = num2str(n);

  fprintf('Computing QNMs for n = %d ... ', n);

  mp.Digits(n);

  m0name = strcat('data/M0_', nstr, '.mat');
  m1name = strcat('data/M1_', nstr, '.mat');
  m2name = strcat('data/M2_', nstr, '.mat');

  M0 = mp.read(m0name);
  M1 = mp.read(m1name);
  M2 = mp.read(m2name);

  e = polyeig(M0, mp('1i')*M1, M2);

  qnms{idx} = e;
  legs{idx} = strcat(' ', nstr);

  fprintf('done.\n');

end % for

%%% Plotting:

figure;

for idx = 1:L
    e = qnms{idx};
    plot(real(e), imag(e), Marks{idx}), grid off, hold on
end % for
hold off

l = legend(legs{:});
set(l, 'box', 'off');
set(l, 'location', 'best');
xlabel('$\mathrm{Re}\,\omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
ylabel('$\mathrm{Im}\,\omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
set(gcf, 'color', 'w');