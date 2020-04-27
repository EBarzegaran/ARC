function [brod_label, lobe_label, hem_label, aal_label, CO] = ParcellateSource(Nodes)
    load Parcellation_precise;
    load Source_CO_3005_nocel;
    load ds3000to400;
    SourceSel = Sources(Nodes);
    CO = CO(SourceSel,:);
    brod_label = brod_label(SourceSel);
    lobe_label = lobe_label(SourceSel);
    hem_label = hem_label(SourceSel);
    load aallabels;
    AAL = [aal_tags.Name aal_tags.Hemi];
    aal_label = AAL(aal_label(SourceSel),:);
end
