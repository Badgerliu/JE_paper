# Single-cell RNA Analysis Data Flow Documentation 

# For Figure 5 scRNA-seq data

# (Arranged by Execution Order)

## 3_Integrate_1st_JE_epi_harmony_imputation_3D.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_20221128.rds`

### Output RDS
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_3D.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE7_20221128.rds`
- `./SeuratObject/JE_combined_epi_MAGIC_GE4_20221128.rds`
- Files in `./scVeloInput/` directory for RNA velocity analysis

---

## 4_Integrate_1st_JE_epi_harmony_imputation_pseudotime_slingshot.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds` (Output from 3_Integrate)

### Output RDS
- `./SeuratObject/JE_combined_epi_harmony_SCT_renamed_slingshot.rds`
- `./SeuratObject/JE_combined_epi_harmony_SCT_renamed_renamed_sce.Rdata`

---

## 7_Integrate_1st_JE_epi_harmony_imputation_pseudotime_slingshot.Rmd

### Input Data
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds` (Output from 3_Integrate)
- `./SeuratObject/JE_combined_epi_harmony_SCT_renamed_slingshot.rds` (Output from 4_Integrate)

### Output RDS
- `./SeuratObject/JE_combined_epi_harmony_SCT_renamed_slingshot_tradeSeq.Rdata`

---

## Figure_5_plots.Rmd

### Input Data
- `./Figure5_scRNA_part/Figure5_scRNA_part.Rdata` (Cache data)
- `./SeuratObject/JE_combined_epi_harmony_SCT_renamed_slingshot_tradeSeq.Rdata` (From 7_Integrate)
- `./SeuratObject/JE_combined_epi_harmony_MAGIC_renamed.rds` (From 3_Integrate)
- `./SeuratObject/JE_combined_epi_harmony_SCT_renamed_renamed_sce.Rdata` (From 4_Integrate)

### Output RDS
- `./Figure5_scRNA_part/Figure5_scRNA_part.Rdata` (Updated cache data)
- Various plot files in `./Figure5_scRNA_part/` directory

---

## Data Flow Diagram
```mermaid
%%{init: {'flowchart': {'nodeSpacing': 100, 'rankSpacing': 150}}}%%
graph TB
    %% Data flow for 3_Integrate
    subgraph sg1["3_Integrate_1st_JE_epi<br/>harmony<br/>imputation_3D.Rmd"]
        harmony[JE_combined_epi_harmony_20221128.rds] --> process1[MAGIC<br/>imputation]
        process1 --> magic[JE_combined_epi_harmony_MAGIC_renamed.rds]
        process1 --> magic3d[JE_combined_epi_MAGIC_3D.rds]
        process1 --> magicGE7[JE_combined_epi_MAGIC_GE7_20221128.rds]
        process1 --> magicGE4[JE_combined_epi_MAGIC_GE4_20221128.rds]
        process1 --> scvelo[scVelo files]
    end

    %% Data flow for 4_Integrate
    subgraph sg2["4_Integrate_1st_JE_epi<br/>harmony<br/>imputation_pseudotime_slingshot.Rmd"]
        magic --> process2[Slingshot<br/>pseudotime]
        process2 --> slingshot[JE_combined_epi_harmony_SCT_renamed_slingshot.rds]
        process2 --> sce[JE_combined_epi_harmony_SCT_renamed_renamed_sce.Rdata]
    end

    %% Data flow for 7_Integrate
    subgraph sg3["7_Integrate_1st_JE_epi<br/>harmony<br/>imputation_pseudotime_slingshot.Rmd"]
        magic --> process3[TradeSeq<br/>analysis]
        slingshot --> process3
        process3 --> tradeseq[JE_combined_epi_harmony_SCT_renamed_slingshot_tradeSeq.Rdata]
    end

    %% Figure_5 input integration
    subgraph sg4["Figure_5_plots.Rmd"]
        magic -->|Imputed data| fig5[Figure 5<br/>visualization]
        sce -->|SCE object| fig5
        tradeseq -->|Gene expression<br/>patterns| fig5
        fig5 -->|Output plots| plots[Figure 5 plots]
    end

    %% Style settings
    classDef process fill:#e1f5fe,stroke:#01579b;
    classDef data fill:#f3e5f5,stroke:#4a148c;
    classDef input fill:#fff3e0,stroke:#e65100;
    classDef output fill:#e8f5e9,stroke:#1b5e20;
    classDef integration fill:#e8f5e9,stroke:#1b5e20,color:#000;
    classDef plotting fill:#fff3e0,stroke:#e65100,color:#000;
    
    class process1,process2,process3 process;
    class harmony,magic,slingshot,tradeseq,magic3d,magicGE7,magicGE4,sce data;
    class scvelo input;
    class fig5,plots output;
    class sg1,sg2,sg3 integration;
    class sg4 plotting;

    %% Annotations
    style harmony fill:#bbdefb,stroke:#1976d2;
    style magic fill:#bbdefb,stroke:#1976d2;
    style slingshot fill:#bbdefb,stroke:#1976d2;
    style sce fill:#bbdefb,stroke:#1976d2;
    style tradeseq fill:#bbdefb,stroke:#1976d2;
    style magic3d fill:#bbdefb,stroke:#1976d2;
    style magicGE7 fill:#bbdefb,stroke:#1976d2;
    style magicGE4 fill:#bbdefb,stroke:#1976d2;
    style plots fill:#e8f5e9,stroke:#1b5e20;
``` 