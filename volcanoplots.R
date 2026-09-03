

##Adipocyte abdomen volcano plot



sig_genes <- rownames(res_adipocyte_face_ordered)[
  !is.na(res_adipocyte_face_ordered$padj) &
    res_adipocyte_face_ordered$padj < 0.05 &
    abs(res_adipocyte_face_ordered$log2FoldChange) >= 1
]

# Volcano plot
EnhancedVolcano(res_adipocyte_face_ordered,
                lab = rownames(res_adipocyte_face_ordered),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                selectLab = sig_genes,
                labSize = 3.5,
                title = "Adipocyte face (keloid vs normal)"
)

##Adipocyte face volcano plot

sig_genes <- rownames(res_keloid_abdomen)[
  !is.na(res_keloid_abdomen$padj) &
    res_keloid_abdomen$padj < 0.05 &
    abs(res_keloid_abdomen$log2FoldChange) >= 1
]

# Volcano plot
EnhancedVolcano(res_keloid_abdomen,
                lab = rownames(res_keloid_abdomen),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                selectLab = sig_genes,
                labSize = 3.5,
                title = "Keloid fibroblasts vs adipocytes abdomen"
)


##Adipocyte fkeloid  volcano plot

sig_genes <- rownames(res_keloid_face)[
  !is.na(res_keloid_face$padj) &
    res_keloid_face$padj < 0.05 &
    abs(res_keloid_face$log2FoldChange) >= 1
]

EnhancedVolcano(res_keloid_face,
                lab = rownames(res_keloid_face),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                selectLab = sig_genes,
                labSize = 3.5,
                title = "keloid fibroblasts vs adipocytes face")




















##Adipocyte normal abdomen volcano plot

sig_genes <- rownames(res_normal_abdomen)[
  !is.na(res_normal_abdomen$padj) &
    res_normal_abdomen$padj < 0.05 &
    abs(res_normal_abdomen$log2FoldChange) >= 1
]

# Volcano plot
EnhancedVolcano(res_normal_abdomen,
                lab = rownames(res_normal_abdomen),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                selectLab = sig_genes,
                labSize = 3.5,
                title = "Normal fibroblasts vs adipocytes abdomen"
)


##Adipocyte normal face volcano plot

sig_genes <- rownames(res_keloid_face)[
  !is.na(res_keloid_face$padj) &
    res_keloid_face$padj < 0.05 &
    abs(res_keloid_face$log2FoldChange) >= 1
]

EnhancedVolcano(res_keloid_face,
                lab = rownames(res_keloid_face),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                selectLab = sig_genes,
                labSize = 3.5,
                title = "keloid fibroblasts vs adipocytes face")

