# Single-cell RNA Analysis Data Flow Documentation 

# For Figure 2 JE scRNA-seq data

# (Arranged by Execution Order)

## 2_Integrate_1st_JE_epi_harmony_subcelltype.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_raw_20221130.rds`

### Output RDS
- Files in `./scVeloInput/` directory for RNA velocity analysis

---

## 3_Integrate_1st_JE_epi_harmony_imputation_3D.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_20221128.rds`

### Output RDS
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_3D.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE7_20221128.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE4_20221128.rds`

---

## 4_Integrate_1st_JE_epi_harmony_imputation_pseudotime_slingshot.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE7_20221128.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE4_20221128.rds`

### Output RDS
- Slingshot pseudotime analysis results and trajectory data

---

## 5_JE_K5_mTmG.Rmd

### Input Data
- `./SeuratObject/JE_K5_D0_Seurat.rds`
- `./SeuratObject/JE_K5_D3_Seurat.rds`
- `./SeuratObject/JE_K5_D5_Seurat.rds`

### Output RDS
- `./SeuratObject/JE_K5_combined.rds`
- `./SeuratObject/JE_K5_combined_harmony.rds`

---

## 6_JE_K5_mTmG_epi.Rmd

### Input Data
- `./SeuratObject/JE_K5_combined_harmony.rds`

### Output RDS
- `./SeuratObject/JE_K5_combined_epi_harmony.rds`
- Epithelial population specific analysis results

---

## Figure_2_plots.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds`
- `./SeuratObject/JE_K5_combined_epi_harmony.rds`
- Slingshot pseudotime analysis results
- RNA velocity analysis results

### Output
- Figure 2 visualization plots and statistical analysis results

---

## Data Flow Diagram
```mermaid
%%{init: {'flowchart': {'nodeSpacing': 100, 'rankSpacing': 150}}}%%
graph TB
    %% Data flow for initial seurat files
    subgraph sg0["JE_K5mTmG Seurat Object Generation"]
        rmd1[JE_K5mTmG_D0_seurat.Rmd] --> d0[JE_K5_D0_Seurat.rds]
        rmd2[JE_K5mTmG_D3_seurat.Rmd] --> d3[JE_K5_D3_Seurat.rds]
        rmd3[JE_K5mTmG_D5_seurat.Rmd] --> d5[JE_K5_D5_Seurat.rds]
    end

    %% Data flow for 2_Integrate
    subgraph sg1["2_Integrate_1st_JE_epi_harmony_subcelltype.Rmd"]
        raw1[JE_combined_epi_raw_20221130.rds] --> process1[Subpopulation<br/>analysis]
        process1 --> scvelo[scVeloInput files]
    end

    %% Data flow for 3_Integrate
    subgraph sg2["3_Integrate_1st_JE_epi_harmony_imputation_3D.Rmd"]
        harmony1[JE_combined_epi_harmony_20221128.rds] --> process2[MAGIC<br/>imputation]
        process2 --> magic1[JE_combined_epi_harmony_MAGIC_renamed.rds]
        process2 --> magic2[JE_combined_epi_MAGIC_GE7_20221128.rds]
        process2 --> magic3[JE_combined_epi_MAGIC_GE4_20221128.rds]
        process2 --> magic4[JE_combined_epi_MAGIC_3D.rds]
    end

    %% Data flow for 4_Integrate
    subgraph sg3["4_Integrate_1st_JE_epi_harmony_imputation_pseudotime_slingshot.Rmd"]
        magic1 & magic2 & magic3 --> process3[Slingshot<br/>pseudotime]
        process3 --> trajectory[Trajectory data]
    end

    %% Data flow for 5_JE_K5_mTmG
    subgraph sg4["5_JE_K5_mTmG.Rmd"]
        d0 & d3 & d5 --> process4[Integration<br/>and harmony]
        process4 --> harmony2[JE_K5_combined.rds]
        process4 --> harmony3[JE_K5_combined_harmony.rds]
        process4 --> epi_raw[JE_K5_epi_raw_20230728.rds]
        magic_harmony[JE_K5_epi_MAGIC_harmony_renamed.rds] --> process4b[Final<br/>processing]
        process4b --> final_harmony[JE_K5_harmony_celltype_20240708.rds]
    end

    %% Data flow for 6_JE_K5_mTmG_epi
    subgraph sg5["6_JE_K5_mTmG_epi.Rmd"]
        epi_raw --> process5[Epithelial<br/>analysis]
        process5 --> magic_harmony
        process5 --> epi[JE_K5_combined_epi_harmony.rds]
    end

    %% Figure_2 input integration
    subgraph sg6["Figure_2_plots.Rmd"]
        magic1 -->|MAGIC data| fig2[Figure 2<br/>visualization]
        epi -->|K5 lineage data| fig2
        trajectory -->|Pseudotime data| fig2
        scvelo -->|Velocity data| fig2
        final_harmony -->|Celltype data| fig2
    end

    %% Style settings
    classDef process fill:#e1f5fe,stroke:#01579b;
    classDef data fill:#f3e5f5,stroke:#4a148c;
    classDef input fill:#fff3e0,stroke:#e65100;
    classDef output fill:#e8f5e9,stroke:#1b5e20;
    
    %% Subgraph styles
    classDef seurat_processing fill:#f5f5f5,stroke:#333333;
    classDef magic_analysis fill:#e8f5e9,stroke:#1b5e20;
    classDef figure_plotting fill:#fff3e0,stroke:#e65100;
    
    class process1,process2,process3,process4,process4b,process5 process;
    class magic1,magic2,magic3,magic4,harmony1,harmony2,harmony3,trajectory,scvelo,epi,epi_raw,magic_harmony,final_harmony data;
    class raw1,d0,d3,d5 input;
    class fig2 output;
    class rmd1,rmd2,rmd3 input;
    
    %% Apply subgraph styles
    class sg0,sg1,sg4,sg5 seurat_processing;
    class sg2,sg3 magic_analysis;
    class sg6 figure_plotting;

    %% Annotations
    style magic1 fill:#bbdefb,stroke:#1976d2;
    style magic2 fill:#bbdefb,stroke:#1976d2;
    style magic3 fill:#bbdefb,stroke:#1976d2;
    style magic4 fill:#bbdefb,stroke:#1976d2;
    style magic_harmony fill:#bbdefb,stroke:#1976d2;