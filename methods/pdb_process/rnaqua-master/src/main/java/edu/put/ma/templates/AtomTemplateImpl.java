package edu.put.ma.templates;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import lombok.Getter;

public class AtomTemplateImpl implements AtomTemplate {

    private List<String> atomNames = null;

    @Getter
    private String templateAtomName = null;

    @Getter
    private String originalAtomName = null;

    public AtomTemplateImpl(final String atomNames, final String templateAtomName) {
        this.originalAtomName = templateAtomName;
        this.templateAtomName = StringUtils.trim(templateAtomName);
        this.atomNames = Arrays.asList(StringUtils.split(atomNames, ','));
    }

    @Override
    public final boolean isProperAtomName(final String patomName) {
        return atomNames.contains(patomName.trim().toUpperCase());
    }

}
