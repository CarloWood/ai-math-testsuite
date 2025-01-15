#include "sys.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Window.h"
#include "debug.h"
#include "sys.h"
#include "utils/AIAlert.h"
#include "utils/ColorPool.h"
#include "utils/debug_ostream_operators.h"
#include <chrono>
#include <random>
#include <sstream>

#if 0
int main(int argc, char *argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  // Handle random seed.
  unsigned int seed;
  if (argc > 1)
  {
    std::istringstream ss(argv[1]);
    ss >> seed;
  }
  else
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  Dout(dc::notice, "Seed: " << seed);

  // Use the mersenne twister engine.
  std::mt19937 engine(seed);

  std::uniform_real_distribution<double> coefficient_dist(-1.0, 1.0);

  // Number of buckets.
  constexpr int N = 1000;

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Epsilon Sums", 1200, 1000);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

    // Create another layer.
    auto layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("layer"));

    // Open the window and start drawing.
    std::thread event_loop([&]() {
      Debug(NAMESPACE_DEBUG::init_thread("event_loop"));
      // Open window, handle event loop. This must be constructed after the draw
      // stuff, so that it is destructed first! Upon destruction it blocks until
      // the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    plot::Plot plot(window.geometry(), {.grid = {.color = color::orange}}, "Probability Distribution", {}, "sum", {},
                    "count", {});
    plot.set_xrange({-1.0, 1.0});
    plot.set_yrange({0.0, 2500.0});
    plot.add_to(background_layer, false);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle dashed_line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});

    // Draw a slider.
    auto slider_r0 = plot.create_slider(layer, {1128, 83, 7, 400}, 0.0, 0.0, 10.0);
    auto slider_r0_label = plot.create_text(layer, slider_style, Pixel{1128, 483}, "r₀");

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in
      // order to avoid flickering.
      window.set_send_expose_events(false);

      // Extract slider values.
      double r0 = slider_r0.value();

      double sum_min = -r0 - 2.0;
      double sum_max = r0 + 2.0;

      std::array<int, N + 1> buckets = {};

      for (int i = 0; i < 1000000; ++i)
      {
        std::array<double, 3> random_numbers;
        for (int j = 0; j < 3; ++j)
          random_numbers[j] = coefficient_dist(engine);
        double sum = random_numbers[0] + r0 * random_numbers[1] + random_numbers[2];
        // Using N buckets.
        int bucket = (sum - sum_min) / (sum_max - sum_min) * N;
        ++buckets[bucket];
      }

      std::array<plot::LinePiece, N + 1> line_pieces;

      for (int b = 0; b <= N; ++b)
      {
        double x = -1.0 + 2.0 / N * b;
        std::cout << x << " : " << buckets[b] << '\n';
        line_pieces[b] = LinePiece{{x, 0.0}, {x, static_cast<double>(buckets[b])}};
        plot.add_line(layer, solid_line_style, plot::LineExtend::none, line_pieces[b]);
      }

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until a redraw is necessary (for example because the user moved a
      // draggable object, or wants to print the current drawing) then go to the
      // top of loop for a redraw.
      if (!window.handle_input_events())
        break; // Program must be terminated.
    }

    event_loop.join();
  } catch (AIAlert::Error const &error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

constexpr size_t n = 4;                 // Dimension of vectors.
constexpr size_t num_samples = 1000000; // Number of Monte Carlo samples.
constexpr size_t num_trials = 100;      // Number of different random vectors a to try.
constexpr double target_length = 1.0;   // Target length for vector a normalization.

// Generate a random unit vector.
std::vector<double> generate_random_vector(std::mt19937& gen)
{
  std::normal_distribution<double> dist(0.0, 1.0);
  std::vector<double> vec(n);

  // Generate random components.
  for (int i = 0; i < n; ++i)
    vec[i] = dist(gen);

  // Normalize to target length.
  double length = std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0));
  for (double& x : vec)
    x *= target_length / length;

  return vec;
}

// Calculate dot product with random uniform vector.
double sample_dot_product(const std::vector<double> &a, std::mt19937 &gen)
{
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  double sum = 0.0;
  for (size_t i = 0; i < n; ++i)
    sum += a[i] * dist(gen);
  return sum;
}

double calculate_variance(std::vector<double> const& arr)
{
  double n = 0;
  double mean = 0;
  double M2 = 0;

  for (double const& value : arr)
  {
    ++n;
    double delta = value - mean;
    mean += delta / n;
    double delta2 = value - mean;
    M2 += delta * delta2;
  }

  return M2 / (n - 1);
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());

  // For storing statistics across trials
  std::vector<double> range_sqrt_var_ratio;

  for (size_t trial = 0; trial < num_trials; ++trial)
  {
    // Generate random normalized vector a.
    auto a = generate_random_vector(gen);       // {a₀, a₁, a₂, a₃, a₄, .. }

#if 0
    std::vector<double> a;
    for (int i = 0; i < n; ++i)
      a.push_back(1.0 / std::sqrt(n));  // a = (1/√n, 1/√n, 1/√n, ...)
#endif

    // Since the length of a is 1, the dot product with the random uniform vector (one sample) will have a variance equal to 1/3.

    // Generate samples.
    std::vector<double> samples;
    samples.reserve(num_samples);
    for (size_t i = 0; i < num_samples; ++i)
      samples.push_back(sample_dot_product(a, gen));

    std::cout << "The variance of the samples is: " << calculate_variance(samples) << '\n';

    // Sort samples to find 25th and 75th percentiles.
    std::sort(samples.begin(), samples.end());
    double p25 = samples[num_samples / 4];
    double p75 = samples[3 * num_samples / 4];
    // Calculate the empirical range containing 50% of the samples.
    double empirical_range = p75 - p25;

    // Store the ratio between the empirical 50% range and the sqrt of the variance (1/3).
    range_sqrt_var_ratio.push_back(empirical_range * std::sqrt(3.0));

    // Assume one coordinate of a dominates, say a = (1, 0, 0, ...).
    // Then the variance of the dot product is ||a||²/3 = 1/3.
    //
    // aₘ²/3 = n aₖ²/3 → |aₘ| = √n |aₖ|
    //
    // while the 50% range ([-b, b]) for both cases become comparible:
    //
    // 1. b ≈ (1² + 0² + 0² + ...)/2 = 0.5
    // and (for large n, where the central limit theorem applies and we add n random variables distributed uniformly in [-1/√n, 1/√n])
    // 2. b ≈ 1.34898 |1/√n| √(n/3) = 1.34898 (1/√n) √(n/3) = 1.34898 / √3 = 0.7788
    //
    // So the ratio between the empirical 50% range and the sqrt of the variance (1/3) should be around 1.5576.

    // Print results for this trial.
    std::cout << "Trial " << trial + 1 << ":\n";
    std::cout << "  Vector a: ";
    for (double x : a)
      std::cout << std::fixed << std::setprecision(4) << x << " ";
    std::cout << "; empirical 50% range: " << empirical_range << "\n";
  }

  // Calculate statistics across all trials.
  double avg_range = std::accumulate(range_sqrt_var_ratio.begin(), range_sqrt_var_ratio.end(), 0.0) / num_trials;
  double var_range = 0.0;
  for (double range : range_sqrt_var_ratio)
    var_range += (range - avg_range) * (range - avg_range);
  var_range /= (num_trials - 1);

  std::cout << "Summary across " << num_trials << " trials:\n";
  std::cout << "Average 50% range: " << avg_range << "\n";
  std::cout << "Standard deviation of range: " << std::sqrt(var_range) << "\n";
}
