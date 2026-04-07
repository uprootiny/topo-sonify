use crate::audio::filter::{FilterMode, TopologicalFilter};

/// SimplicialPluck: a Karplus-Strong string synthesizer where the
/// string's overtone structure is shaped by the simplicial complex.
///
/// ## Mathematical Foundation
///
/// The Karplus-Strong algorithm models a vibrating string as a delay line
/// with a lowpass filter in the feedback loop. This is a discrete
/// approximation to the 1D wave equation ∂²u/∂t² = c² ∂²u/∂x².
///
/// The eigenmodes of the wave equation on an interval [0, L] are
/// uₙ(x,t) = sin(nπx/L) · cos(nπct/L), with frequencies fₙ = nc/(2L).
///
/// On a simplicial complex K, we replace the 1D Laplacian with the
/// combinatorial Laplacian Lₖ = ∂ₖ₊₁ ∂ₖ₊₁ᵀ + ∂ₖᵀ ∂ₖ.
/// The eigenvalues of Lₖ determine the "natural frequencies" of the
/// k-dimensional "string" on K.
///
/// ## Design
///
/// - The delay line length = sample_rate / frequency (standard KS)
/// - The f-vector shapes the initial excitation waveform:
///   instead of white noise, we use a signal whose spectrum
///   has energy concentrated at f-vector-derived harmonics.
/// - The Betti numbers control the feedback filter:
///   β₁ > 0 → more resonant feedback (loops sustain energy)
///   β₀ > 1 → split the delay into parallel lines (detuned)
pub struct SimplicialPluck {
    delay_lines: Vec<Vec<f64>>,
    positions: Vec<usize>,
    filters: Vec<TopologicalFilter>,
    sample_rate: f64,
    damping: f64,
    active: bool,
}

impl SimplicialPluck {
    pub fn new(frequency: f64, sample_rate: f64) -> Self {
        let delay_len = (sample_rate / frequency) as usize;
        Self {
            delay_lines: vec![vec![0.0; delay_len]],
            positions: vec![0],
            filters: vec![TopologicalFilter::new(
                FilterMode::LowPass,
                frequency * 4.0,
                0.707,
                sample_rate,
            )],
            sample_rate,
            damping: 0.996,
            active: false,
        }
    }

    /// Configure from a simplicial complex's properties.
    pub fn configure(&mut self, frequency: f64, f_vector: &[usize], b0: usize, b1: usize) {
        let num_strings = b0.clamp(1, 4);
        let delay_len = (self.sample_rate / frequency) as usize;

        self.delay_lines.clear();
        self.positions.clear();
        self.filters.clear();

        for i in 0..num_strings {
            let detune = 1.0 + (i as f64 * 0.003);
            let len = (delay_len as f64 / detune) as usize;
            self.delay_lines.push(vec![0.0; len.max(2)]);
            self.positions.push(0);

            // Higher β₁ → higher resonance in the feedback filter
            let q = 0.5 + b1 as f64 * 0.3;
            self.filters.push(TopologicalFilter::new(
                FilterMode::LowPass,
                frequency * (2.0 + f_vector.len() as f64),
                q.min(5.0),
                self.sample_rate,
            ));
        }

        // Damping decreases (more sustain) with more loops
        self.damping = (0.999 - b1 as f64 * 0.0005).clamp(0.98, 0.9999);
    }

    /// Excite the string — "pluck" it.
    /// Uses f-vector-weighted noise instead of pure white noise.
    pub fn pluck(&mut self, f_vector: &[usize]) {
        self.active = true;
        let total: usize = f_vector.iter().sum::<usize>().max(1);

        for line in &mut self.delay_lines {
            let line_len = line.len();
            for (i, sample) in line.iter_mut().enumerate() {
                // Excitation signal: weighted sum of harmonics from f-vector
                let mut val = 0.0;
                for (k, &fk) in f_vector.iter().enumerate() {
                    let weight = fk as f64 / total as f64;
                    let harmonic = (k + 1) as f64;
                    val += weight
                        * (std::f64::consts::TAU * harmonic * i as f64 / line_len as f64).sin();
                }
                // Add a touch of noise for naturalness
                let noise = ((i * 1103515245 + 12345) % 65536) as f64 / 32768.0 - 1.0;
                *sample = val * 0.7 + noise * 0.3;
            }
        }
    }

    pub fn next_sample(&mut self) -> f64 {
        if !self.active {
            return 0.0;
        }

        let mut output = 0.0;
        let num = self.delay_lines.len();

        for i in 0..num {
            let pos = self.positions[i];
            let len = self.delay_lines[i].len();
            let next_pos = (pos + 1) % len;

            // Karplus-Strong: average adjacent samples and apply damping
            let current = self.delay_lines[i][pos];
            let next = self.delay_lines[i][next_pos];
            let averaged = (current + next) * 0.5 * self.damping;

            // Apply feedback filter
            let filtered = self.filters[i].process(averaged);
            self.delay_lines[i][pos] = filtered;

            output += current;
            self.positions[i] = next_pos;
        }

        output / num as f64
    }

    pub fn fill_buffer(&mut self, buffer: &mut [f64]) {
        for sample in buffer.iter_mut() {
            *sample = self.next_sample();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pluck_produces_sound() {
        let mut pluck = SimplicialPluck::new(440.0, 44100.0);
        pluck.pluck(&[4, 6, 4]); // tetrahedron f-vector
        let mut buffer = vec![0.0; 1024];
        pluck.fill_buffer(&mut buffer);
        let energy: f64 = buffer.iter().map(|s| s * s).sum();
        assert!(energy > 0.0);
    }
}
