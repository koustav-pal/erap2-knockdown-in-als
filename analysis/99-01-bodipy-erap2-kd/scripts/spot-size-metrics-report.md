# Report: BODIPY spot-size metrics in the ERAP2 knockdown analysis

## Scope

This report explains how the following metrics are calculated in `03-plot-bodipy-spot-area.R`:

1. Relative spot size
2. Average spot size scaled to 100 µm² of cytoplasm

The analysis combines individual-spot measurements from the spot files with total spot area and cytoplasmic area measurements from the corresponding cell-information files.

## 1. Spot area scaled to 100 µm²

For each cell, the total spot area is normalised by the cell's cytoplasmic area:

\[
\text{total spot area scaled to 100 µm²}
=
\frac{\text{total spot area in the cell (µm²)}}
{\text{cytoplasmic area of the cell (µm²)}}
\times 100
\]

This measures how much spot area would be expected in a standardised 100 µm² region of cytoplasm if the observed spot density remained constant.

Interpretation:

- A value of 5 means that spots occupy an area equivalent to 5 µm² per 100 µm² of cytoplasm.
- This is a cytoplasm-normalised measure of total spot coverage, not the physical area of one spot.
- The script and plot labels use “100 um,” but the correct area unit is **100 µm²**.

Relevant code:

```r
scaled_spot_area <- (spot_area_in_um / cytoplasmic_area) * 100
```

## 2. Average spot size scaled to 100 µm²

The script calculates a per-cell mean by dividing the normalised total spot area by the number of spots in that cell:

\[
\text{mean scaled spot size}
=
\frac{\text{total spot area scaled to 100 µm²}}
{\text{number of spots in the cell}}
\]

Equivalently:

\[
\text{mean scaled spot size}
=
\frac{\text{total spot area (µm²)}}
{\text{cytoplasmic area (µm²)}}
\times
\frac{100}{\text{number of spots}}
\]

Relevant code:

```r
summarised_cell_df$mean_size_of_a_spot <-
    summarised_cell_df$spot_area_scaled_to_100um /
    summarised_cell_df$number_of_spots
```

Interpretation:

- This is the mean spot area per cell after normalising for cytoplasmic area.
- It is useful for comparing cells of different sizes.
- It is not the unnormalised physical mean spot area. That would instead be `total spot area / number of spots`.

At the cell-line/treatment level, the script reports the **median** of these per-cell means. Therefore, plot labels using “average” can be misleading: the spot size is first calculated as a mean within each cell and then summarised using the median across cells.

## 3. Individual spot area

The cell-information file provides total spot area in µm², while the individual spot file provides each spot's area in pixels. The script estimates the physical area of each individual spot by allocating the cell's total measured spot area according to the spot's proportion of total pixel area:

\[
\text{individual spot area (µm²)}
=
\frac{\text{individual spot area (pixels)}}
{\text{sum of spot areas in that cell (pixels)}}
\times
\text{total spot area in the cell (µm²)}
\]

The individual value is then scaled to 100 µm² of cytoplasm:

\[
\text{individual scaled spot area}
=
\frac{\text{individual spot area (µm²)}}
{\text{cytoplasmic area (µm²)}}
\times 100
\]

## 4. Relative spot size

For each individual spot, relative spot size is calculated by dividing its scaled area by the mean scaled spot size of all spots in the same cell:

\[
\text{relative spot size}
=
\frac{\text{individual scaled spot area}}
{\text{mean scaled spot size in the same cell}}
\]

Because both numerator and denominator contain the same cytoplasmic-area scaling factor, this simplifies to:

\[
\text{relative spot size}
=
\frac{\text{individual spot area}}
{\text{mean spot area in the same cell}}
\]

Relevant code:

```r
bodypy_df$relative_spot_size <-
    bodypy_df$spot_area_scaled_to_100um /
    summarised_cell_df$mean_size_of_a_spot[
        match(bodypy_df$temp_id, summarised_cell_df$temp_id)
    ]
```

Interpretation:

- Relative spot size is dimensionless.
- `1` means the spot is the mean size for its cell.
- `2` means the spot is twice the mean size for its cell.
- `0.5` means the spot is half the mean size for its cell.
- The mean relative spot size within a cell is necessarily 1, subject to numerical precision.

For downstream cell-level analysis, the script retains the **maximum relative spot size in each cell**:

```r
relative_spot_size = max(x$relative_spot_size)
```

Thus, the downstream `relative_spot_size` variable represents the size of the largest spot relative to that cell's mean spot size. It does not represent the average relative size of all spots.

## 5. Sampling strategy

The downstream analysis uses repeated, probability-weighted sampling to account for differences in cytoplasmic area among experimental fields. Sampling is performed after retaining the `ERAP2`, `Scramble`, and `Untx` treatments and excluding the `CTRL3` cell line.

### 5.1 Field-level cytoplasmic-area weights

Each experimental field is identified by the combination of experiment file, row, column, and time point. The total cytoplasmic area of every field is divided by the total cytoplasmic area across all retained fields:

\[
P_{\mathrm{field}}
=
\frac{\text{total cytoplasmic area of the field}}
{\text{total cytoplasmic area across all fields}}
\]

This value is stored as `sample_total_cyto_area_probability` for every cell belonging to that field.

The script also calculates `cell_selection_probability` using individual-cell cytoplasmic area, total cytoplasmic area, and the original field sampling proportion. However, this variable is **not used** by the subsequent `sample()` call. The implemented sampling weights are therefore field-level total-cytoplasmic-area weights repeated for each cell, rather than the calculated individual-cell selection probabilities.

### 5.2 Repeated sampling

Sampling is carried out as follows:

1. Split the complete cell dataset by treatment.
2. For each treatment, sample 1,000 cell-row indices using `sample_total_cyto_area_probability` as the probability weights.
3. Because `replace` is not specified in the R `sample()` call, sampling is performed **without replacement**.
4. Calculate the proportion of selected cells originating from each experimental field.
5. Compare the ranks of the selected field proportions with the ranks of the original field-level cell proportions using a Wilcoxon test.
6. Reject the sample and draw another while the Wilcoxon p-value is below 0.05. Accept the sample once p ≥ 0.05.
7. Repeat this procedure to obtain 1,000 accepted samples for every treatment.

Relevant code:

```r
while (pval < 0.05) {
    current_vector <- sample(
        x = seq_along(treatment_df$sample_total_cyto_area_probability),
        size = 1000,
        prob = treatment_df$sample_total_cyto_area_probability
    )

    current_treatment_df <- treatment_df[current_vector, ]
    current_proportions <- table(
        paste(
            current_treatment_df$experiment,
            current_treatment_df$experiment_label,
            sep = "_"
        )
    ) / nrow(current_treatment_df)

    pval <- wilcox.test(
        rank(current_proportions),
        rank(sample_cell_proportion[names(current_proportions)])
    )$p.value
}
```

### 5.3 Statistics calculated from each accepted sample

Within every accepted treatment-specific sample, cells are divided into `control` and `vcp` conditions. The following medians are calculated for each condition:

- `median_spot_size`: median of the per-cell mean spot size scaled to 100 µm²;
- `spot_area_scaled_to_100um`: median of the per-cell median individual spot area scaled to 100 µm²;
- `median_relative_spot_size`: median of the largest-spot-to-mean-spot ratio from each cell.

This produces one summary for each treatment-condition combination in each of the 1,000 iterations. The distributions of these 1,000 summaries are then shown using boxplots and quasirandom points. Wilcoxon tests are applied directly to the iteration-level summaries to compare conditions within treatments and treatments within conditions.

### 5.4 Methodological cautions

- The code requires at least 1,000 available cell rows in every treatment because sampling is performed without replacement.
- No maximum number of attempts is specified for the `while` loop. Sampling could therefore continue indefinitely if the acceptance criterion cannot be met.
- No random seed is set in the script, so the sampled results are not exactly reproducible between runs.
- The field-level weight is repeated for every cell in a field. Consequently, a field's effective total sampling weight depends on both its total cytoplasmic area and its number of eligible cells.
- The Wilcoxon acceptance test checks whether the ranked field proportions differ detectably; it does not demonstrate that the selected and original field proportions are equivalent.
- The 1,000 iteration-level summaries come from repeated samples of the same underlying cells and should not automatically be interpreted as 1,000 independent biological replicates.

## 6. Important naming issue in the downstream summary

The `max_spot_size_per_cell` table contains:

```r
spot_area_scaled_to_100um = median(x$spot_area_scaled_to_100um)
relative_spot_size = max(x$relative_spot_size)
```

Consequently:

- `spot_area_scaled_to_100um` in this table is the **median individual scaled spot area per cell**.
- `relative_spot_size` is the **maximum relative spot size per cell**.

These two columns summarise different parts of the within-cell distribution and should not both be described simply as “average spot size.”

## 7. Recommended terminology

To make the analysis and figure labels clearer, use:

| Current variable or label | Recommended description |
|---|---|
| `scaled_spot_area` at cell level | Total spot area per 100 µm² of cytoplasm |
| `mean_size_of_a_spot` | Mean spot area per spot, normalised to 100 µm² of cytoplasm |
| Individual `spot_area_scaled_to_100um` | Individual spot area normalised to 100 µm² of cytoplasm |
| Downstream `spot_area_scaled_to_100um` | Median individual spot area per cell, normalised to 100 µm² |
| Individual `relative_spot_size` | Spot area relative to the mean spot area in the same cell |
| Downstream `relative_spot_size` | Largest-spot-to-mean-spot area ratio per cell |

## 8. Summary

Relative spot size describes how large an individual spot is compared with the mean spot in the same cell. The downstream analysis specifically uses the largest such ratio per cell. Average spot size scaled to 100 µm² describes the mean spot area after accounting for cytoplasmic area, enabling comparisons across differently sized cells. Group-level figures use medians across sampled cells rather than arithmetic means. The sampling procedure generates 1,000 accepted, treatment-specific samples of 1,000 cells using field-level cytoplasmic-area weights and a Wilcoxon-based field-composition acceptance rule.
