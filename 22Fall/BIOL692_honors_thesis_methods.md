# Identification of Repressive RNAs with Xist-Like Functions in the Mouse Transcriptome: Materials and Methods

## Quantification of RNA-Protein Interactions
Mickey Murvin - another member of the Calabrese Lab - performed RNA Immuprecipitation (RIP) on mouse trophoblast stem cells (TSCs)  with antibodies of 27 RNA-binding proteins important for Xist function: `Aly/Ref,G9a,HnrnpC,HNRNPK,HnrnpM,HnrnpU,Jarid2,LBR,MAtr3,Nudt21,PABPN1,PTBP1,RBM15,Ring1b,RYBP,SAFB,SPEN,SRSF1,SUPT16H,SUZ12,U2AF35,XRN2,Tia1,Ciz1,U2AF65,Nxf1,SFPQ`. `IgG` antibodies were used as control.

I stored her RIP data in `/proj/calabrlb/users/Zhiyue/22_sp/rip/` on Longleaf:
```
tsc_b12_alyref_rip_S6_R1_001.fastq.gz
tsc_b12_g9a_rip_S6_R1_001.fastq.gz
tsc_b12_hnrnpc_rip_S4_R1_001.fastq.gz
tsc_b12_hnrnpk_rip_S7_R1_001.fastq.gz
tsc_b12_hnrnpm_rip_S5_R1_001.fastq.gz
tsc_b12_hnrnpu_rip_S3_R1_001.fastq.gz
tsc_b12_igg_goat_rip_S2_R1_001.fastq.gz
tsc_b12_igg_mouse_rip_S1_R1_001.fastq.gz
tsc_b12_igg_rb_rip_S2_R1_001.fastq.gz
tsc_b12_jarid2_rip_S5_R1_001.fastq.gz
tsc_b12_lbr_rip_S7_R1_001.fastq.gz
tsc_b12_matrin3_rip_S3_R1_001.fastq.gz
tsc_b12_pabpn1_rip_S8_R1_001.fastq.gz
tsc_b12_ptbp1_rip_S9_R1_001.fastq.gz
tsc_b12_rbm15_rip_S8_R1_001.fastq.gz
tsc_b12_rybp_rip_S12_R1_001.fastq.gz
tsc_b12_srsf1_rip_S4_R1_001.fastq.gz
tsc_b15_nudt21_b1_rip_S8_R1_001.fastq.gz
tsc_b15_supt16h_b1_rip_S7_R1_001.fastq.gz
tsc_b15_u2af35_b2_rip_S13_R1_001.fastq.gz
tsc_b15_u2af65_b2_rip_S16_R1_001.fastq.gz
tsc_b15_xrn2_b2_rip_S14_R1_001.fastq.gz
tsc_B5_alyref-rip_S8_R1_001.fastq.gz
tsc_B5_hnrnpk-rip_S10_R1_001.fastq.gz
tsc_B5_igg_mouse-rip_S2_R1_001.fastq.gz
tsc_B5_matrin3-rip_S14_R1_001.fastq.gz
tsc_B5_rybp-rip_S15_R1_001.fastq.gz
tsc_B5_spen-rip_S13_R1_001.fastq.gz
tsc_b8_igg_goat_rip_S3_R1_001.fastq.gz
tsc_b8_igg_rb_rip_S2_R1_001.fastq.gz
tsc_b8_jarid2_rip_S10_R1_001.fastq.gz
tsc_b8_ring1b_rip_S12_R1_001.fastq.gz
tsc_b8_suz12_rip_S16_R1_001.fastq.gz
tsc_c3_g9a_rip_S2_R1_001.fastq.gz
tsc_c3_hnrnpc_rip_S10_R1_001.fastq.gz
tsc_c3_hnrnpm_rip_S9_R1_001.fastq.gz
tsc_c3_hnrnpu_rip_S8_R1_001.fastq.gz
tsc_c3_lbr_rip_S4_R1_001.fastq.gz
tsc_c3_nudt21_rip_S14_R1_001.fastq.gz
tsc_c3_pabpn1_rip_S3_R1_001.fastq.gz
tsc_c3_ptbp1_rip_S5_R1_001.fastq.gz
tsc_c3_rbm15_rip_S7_R1_001.fastq.gz
tsc_c3_srsf1_rip_S6_R1_001.fastq.gz
tsc_c3_supt16h_rip_S13_R1_001.fastq.gz
tsc_c3_tia1_1_rip_S15_R1_001.fastq.gz
tsc_c3_tia1_2_rip_S16_R1_001.fastq.gz
tsc_c3_u2af35_rip_S12_R1_001.fastq.gz
tsc_c3_xrn2_rip_S11_R1_001.fastq.gz
tsc_c9_ciz1_2_rip_S18_R1_001.fastq.gz
tsc_c9_hnrnpc_rip_S15_R1_001.fastq.gz
tsc_c9_hnrnpm_rip_S14_R1_001.fastq.gz
tsc_c9_hnrnpu_rip_S16_R1_001.fastq.gz
tsc_c9_spen_rip_S21_R1_001.fastq.gz
tsc_c9_srsf1_rip_S19_R1_001.fastq.gz
tsc_c9_u2af35_rip_S20_R1_001.fastq.gz
tsc_nxf1_S8_R1_001.fastq.gz
tsc_ring1b_S7_R1_001.fastq.gz
tsc_safb1-rip_mm_S19_R1_001.fastq
tsc_sfpq_S10_R1_001.fastq.gz
tsc_spen_novus_S15_R1_001.fastq.gz
tsc_suz12_S6_R1_001.fastq.gz
```

