package edu.put.ma.model.services;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "deformation-profile")
public class DeformationProfile {
    @XmlElement(name = "helices")
    private List<Helix> helixList = new ArrayList<Helix>();

    @XmlElement(name = "loops")
    private List<Loop> loopList = new ArrayList<Loop>();

    @XmlElement(name = "draw")
    private List<DrawComponent> drawList = new ArrayList<DrawComponent>();

    @XmlElement
    private Boolean isDataGenerated = false;

    @XmlElement
    private Boolean isImageGenerated = true;

    @XmlElement
    private Boolean isDetectedAutomatically = false;

    public final List<Helix> getHelicesList() {
        return helixList;
    }

    public final void setHelixList(final List<Helix> phelixList) {
        this.helixList = phelixList;
    }

    public final void addHelix(final Helix helix) {
        this.helixList.add(helix);
    }

    public final void addAllHelices(final List<Helix> phelixList) {
        this.helixList.addAll(phelixList);
    }

    public final List<Loop> getLoopsList() {
        return loopList;
    }

    public final void setLoopList(final List<Loop> ploopList) {
        this.loopList = ploopList;
    }

    public final void addLoop(final Loop loop) {
        this.loopList.add(loop);
    }

    public final void addAllLoops(final List<Loop> ploopList) {
        this.loopList.addAll(ploopList);
    }

    public final List<DrawComponent> getDrawsList() {
        return drawList;
    }

    public final void setDrawList(final List<DrawComponent> pdrawList) {
        this.drawList = pdrawList;
    }

    public final void addDraw(final DrawComponent component) {
        this.drawList.add(component);
    }

    public final void addAllDraws(final List<DrawComponent> pdrawList) {
        this.drawList.addAll(pdrawList);
    }

    public final Boolean getIsDataGenerated() {
        return isDataGenerated;
    }

    public final void setIsDataGenerated(final boolean pisDataGenerated) {
        this.isDataGenerated = pisDataGenerated;
    }

    public final Boolean getIsImageGenerated() {
        return isImageGenerated;
    }

    public final void setIsImageGenerated(final boolean pisImageGenerated) {
        this.isImageGenerated = pisImageGenerated;
    }

    public final Boolean getIsDetectedAutomatically() {
        return isDetectedAutomatically;
    }

    public final void setIsDetectedAutomatically(final boolean pisDetectedAutomatically) {
        this.isDetectedAutomatically = pisDetectedAutomatically;
    }
}
