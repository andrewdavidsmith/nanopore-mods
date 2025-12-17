/* MIT License
 *
 * Copyright (c) 2025 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include "CLI11.hpp"
#include "json.hpp"

#include <htslib/sam.h>

#include <array>
#include <cstdint>
#include <fstream>
#include <print>
#include <stdexcept>
#include <string>

// clang-format off
static constexpr std::array<std::uint8_t, 256> encoding = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 16
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 32
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 48
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 64
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 80
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 96
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 112
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 128
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 144
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 160
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 176
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 192
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 208
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 224
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 240
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4   // 256
};
// clang-format on

static constexpr auto n_nucs = 4u;

struct mod_prob_stats {
  static constexpr auto max_mods = 10;
  static constexpr auto n_values = 256;
  // scratch
  std::array<hts_base_mod, max_mods> mods{};
  std::unique_ptr<hts_base_mod_state, void (*)(hts_base_mod_state *)> m;

  std::array<std::array<std::uint64_t, n_values>, n_nucs> methyl_fwd{};
  std::array<std::array<std::uint64_t, n_values>, n_nucs> methyl_rev{};
  std::array<std::array<std::uint64_t, n_values>, n_nucs> hydroxy_fwd{};
  std::array<std::array<std::uint64_t, n_values>, n_nucs> hydroxy_rev{};

  mod_prob_stats() : m{hts_base_mod_state_alloc(), &hts_base_mod_state_free} {};

  [[nodiscard]] auto
  operator()(const bam1_t *aln) {
    static constexpr auto h_idx = 0;
    static constexpr auto m_idx = 1;
    const auto qlen = aln->core.l_qseq;
    const auto seq = bam_get_seq(aln);
    const auto d = mods.data();
    const auto is_rev = bam_is_rev(aln);

    bam_parse_basemod(aln, m.get());
    // ADS: or bam_parse_basemod2(aln, m.get(), HTS_MOD_REPORT_UNCHECKED)

    int pos{};
    int n{};
    while ((n = bam_next_basemod(aln, m.get(), d, max_mods, &pos)) > 0) {
      if (n < m_idx)
        continue;
      const auto other_nuc =
        is_rev ? (pos > 0 ? seq_nt16_str[bam_seqi(seq, pos - 1)] : '\0')
               : (pos + 1 < qlen ? seq_nt16_str[bam_seqi(seq, pos + 1)] : '\0');
      // NOLINTBEGIN(*-constant-array-index)
      const auto other_enc = encoding[static_cast<std::uint8_t>(other_nuc)];
      if (other_enc == n_nucs)
        continue;
      if (is_rev) {
        hydroxy_rev[other_enc][mods[h_idx].qual]++;
        methyl_rev[other_enc][mods[m_idx].qual]++;
      }
      else {
        hydroxy_fwd[other_enc][mods[h_idx].qual]++;
        methyl_fwd[other_enc][mods[m_idx].qual]++;
      }
      // NOLINTEND(*-constant-array-index)
    }
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(mod_prob_stats, methyl_fwd, methyl_rev,
                                 hydroxy_fwd, hydroxy_rev)
};

int
main(int argc, char *argv[]) {  // NOLINT(*-c-arrays)
  std::string outfile;
  std::string infile;

  CLI::App app{};
  argv = app.ensure_utf8(argv);
  app.usage("Usage: nanopore-mods [options]");

  // clang-format off
  app.add_option("-i,--input", infile, "BAM/SAM input file")
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output", outfile, "JSON output file")
    ->required();
  // clang-format on

  if (argc < 2) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  auto in = hts_open(infile.data(), "r");
  if (!in)
    throw std::runtime_error("failed to open file: " + infile);
  std::unique_ptr<sam_hdr_t, void (*)(sam_hdr_t *)> hdr{sam_hdr_read(in),
                                                        &bam_hdr_destroy};
  if (!hdr)
    throw std::runtime_error("failed to parse header from file: " + infile);

  std::unique_ptr<bam1_t, void (*)(bam1_t *)> aln{bam_init1(), &bam_destroy1};

  mod_prob_stats mps;

  std::int32_t read_status{};
  while ((read_status = sam_read1(in, hdr.get(), aln.get())) > -1)
    mps(aln.get());

  hts_close(in);

  if (read_status < -1) {  // -1 is EOF
    std::println(std::cerr, "failed reading bam record");
    return EXIT_FAILURE;
  }

  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("Error opening output file: " + outfile);
  out << nlohmann::json(mps).dump(4) << "\n";

  return EXIT_SUCCESS;
}
