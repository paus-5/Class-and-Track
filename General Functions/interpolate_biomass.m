function OTU_abundance = interpolate_biomass(t_obs,biomass,t_OTU,OTU_rel)
biomass_interpolation = interp1(t_obs,biomass,t_OTU);
OTU_abundance = diag(biomass_interpolation)*OTU_rel;
end