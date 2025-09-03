package edu.put.ma.model.services;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "structures")
public class StructuresSet {
    @XmlElement(name = "atomsFiltering")
    private Boolean isFiltering = false;

    @XmlElement
    private String annotationTool = "MC-Annotate";

    @XmlElement(name = "structure")
    private List<PdbStructure> structureList = new ArrayList<PdbStructure>();

    @XmlElement(name = "deformation-profile")
    private DeformationProfile deformationProfile = new DeformationProfile();

    public final DeformationProfile getDeformationsProfile() {
        return deformationProfile;
    }

    public final void setDeformationProfile(final DeformationProfile pdeformationProfile) {
        this.deformationProfile = pdeformationProfile;
    }

    public final String getAnnotationsTool() {
        return annotationTool;
    }

    public final void setAnnotationTool(final String pannotationTool) {
        this.annotationTool = pannotationTool;
    }

    public final Boolean getIsFilterings() {
        return isFiltering;
    }

    public final void setIsFiltering(final Boolean pisFiltering) {
        this.isFiltering = pisFiltering;
    }

    public final List<PdbStructure> getStructuresList() {
        return structureList;
    }

    public final void setStructureList(final List<PdbStructure> pstructureList) {
        this.structureList = pstructureList;
    }

    public final void addStructure(final PdbStructure str) {
        this.structureList.add(str);
    }

    public final void addStructure(final PdbStructure str, final int index) {
        this.structureList.add(index, str);
    }

    public final void clearStructureList() {
        this.structureList.clear();
    }

    public final PdbStructure getStructureByIndex(final int idx) {
        return this.structureList.get(idx);
    }

    public final int getStructureListSize() {
        return this.structureList.size();
    }

    public final PdbStructure getStructureByFilename(final String filename) {
        for (int i = 0; i < this.structureList.size(); i++) {
            if (this.structureList.get(i).getFilename().equals(filename)) {
                return this.structureList.get(i);
            }
        }
        return null;
    }
}
