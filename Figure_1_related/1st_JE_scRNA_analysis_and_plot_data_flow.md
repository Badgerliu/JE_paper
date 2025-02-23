# Single-cell RNA Analysis Data Flow Documentation 

# For 1st round JE scRNA-seq data

# (Arranged by Execution Order)

## 1_Integrate_1st_JE_harmony_celltype_annotation.Rmd

### Input Data
- `./SeuratObject/JE_D0_raw.rds`
- `./SeuratObject/JE_D3_raw.rds`
- `./SeuratObject/JE_D5_raw.rds`

### Output RDS
- `./SeuratObject/JE_D0_filter.rds`
- `./SeuratObject/JE_D3_filter.rds`
- `./SeuratObject/JE_D5_filter.rds`
- `./SeuratObject/JE_combined.rds`
- `./SeuratObject/JE_combined_harmony.rds`

---

## 2_Integrate_1st_JE_epi_harmony_subcelltype.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_20221128.rds` (Output from 1_Integrate)

### Output RDS
- Intermediate files in `./scVeloInput/` directory (No direct RDS output)

---

## 3_Integrate_1st_JE_epi_harmony_imputation_3D.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_20221128.rds` (Output from 1_Integrate)
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds` (Self-iteratively generated)

### Output RDS
- `./SeuratObject/JE_combined_epi_MAGIC_3D.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE7_20221128.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE4_20221128.rds`

---

## Figure_1_plots.Rmd

### Input Data
- `./SeuratObject/JE_combined_harmony_20221128.rds` (From 1_Integrate)
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds` (From 3_Integrate)
- `./SeuratObject/perio_seurat_20221008.rds` (External reference data)
- `./Figure1_scRNA_part/Figure1_scRNA_part.Rdata` (Self cache)

### Output RDS
- `./Figure1_scRNA_part/Figure1_scRNA_part.Rdata`

---

## Data Flow Diagram
```mermaid
graph TD
    %% Data flow for 1_Integrate
    subgraph "1_Integrate_1st_JE_harmony_celltype_annotation.Rmd"
        raw1[JE_D0_raw.rds] & raw2[JE_D3_raw.rds] & raw3[JE_D5_raw.rds] --> process1[Data filtering and integration]
        process1 --> harmony[JE_combined_harmony_20221128.rds]
    end

    %% Data flow for 2_Integrate
    subgraph "2_Integrate_1st_JE_epi_harmony_subcelltype.Rmd"
        harmony --> process2[Subpopulation analysis]
        process2 --> scvelo[scVelo intermediate files]
    end

    %% Data flow for 3_Integrate
    subgraph "3_Integrate_1st_JE_epi_harmony_imputation_3D.Rmd"
        harmony --> process3[MAGIC imputation]
        process3 --> magic[JE_combined_epi_harmony_MAGIC_renamed.rds]
    end

    %% Figure_1 input integration
    subgraph "Figure_1_plots.Rmd"
        harmony -->|Main integrated data| fig1[Figure_1 visualization]
        magic -->|Imputed data| fig1
        ext[perio_seurat_20221008.rds] -->|External reference data| fig1
        cache[Figure1_scRNA_part.Rdata] -->|Cached data| fig1
    end

    %% Style settings
    classDef process fill:#e1f5fe,stroke:#01579b;
    classDef data fill:#f3e5f5,stroke:#4a148c;
    classDef input fill:#fff3e0,stroke:#e65100;
    classDef output fill:#e8f5e9,stroke:#1b5e20;
    
    class process1,process2,process3 process;
    class harmony,magic,scvelo data;
    class raw1,raw2,raw3,ext input;
    class fig1 output;

    %% Annotations
    style harmony fill:#bbdefb,stroke:#1976d2;
    style magic fill:#bbdefb,stroke:#1976d2;
