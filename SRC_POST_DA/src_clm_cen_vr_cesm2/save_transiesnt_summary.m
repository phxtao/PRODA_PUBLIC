function save_transiesnt_summary(save_pathway, online_tot_eco_c, online_tot_veg_c, online_tot_lit_c, online_npp, online_soc_stock, offline_soc_stock);

transiesnt_summary.online_tot_eco_c = online_tot_eco_c;
transient_summary.online_tot_veg_c = online_tot_veg_c;
transient_summary.online_tot_lit_c = online_tot_lit_c;
transient_summary.online_npp = online_npp;
transient_summary.online_soc_stock = online_soc_stock;
transient_summary.offline_soc_stock = offline_soc_stock;

save(save_pathway, 'transient_summary');
end