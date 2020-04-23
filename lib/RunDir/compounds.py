# python3
# Python Translation of feba/lib/Compounds.pm
import os
import re
import logging
from feba_utils import read_table, read_column_names

# LoadCompounds happens before LoadMedia

# First run
def LoadCompounds(compounds_dir, compounds_dict, synonyms_dict):
    # Uses subfunction SynToKey, read_table, read_column_names 

    compounds_file = os.path.join(compounds_dir, "Compounds.tsv")
    req = "Compound CAS FW Synonyms".split(" ")
    compounds = read_table(compounds_file, req)
    headers = read_column_names(compounds_file)
    if not len(compounds) > 0:
        raise Exception("No rows in " + compounds_file)

    for row in compounds:
        compound = row["Compound"]
        if compound in compounds_dict:
            raise Exception("Duplicate compound id {}".format(compound))
        mw = row["FW"]
        if mw == "NA":
            mw = ""
        if not (mw == "" or re.search(r"^[0-9]+[.]?[0-9]*$", mw)):
            raise Exception(
                "Invalid weight {} for compound {} in {}\n \
                        ".format(
                    mw, compound, compounds_file
                )
            )
        cas = row["CAS"]
        cas.strip()
        if cas == "NA":
            cas = ""
        # Below accounts for cases where cas is multiple CAS divided by '/'
        cas = cas.split("/")[0].strip()
        if not (cas == "" or (re.search(r"^\d+[0-9-]*\d$", cas))):
            raise Exception(
                "Invalid cas number '{}' for compound {} in {}\n \
                        ".format(
                    mw, compound, compounds_file
                )
            )
        compounds_dict[compound] = [compound, cas, mw]

        # Synonyms
        syns = [compound]
        if re.search(r";", row["Synonyms"]):
            # allow semicolon separators instead of comma separators
            syns += row["Synonyms"].split("; ")
        else:
            syns += row["Synonyms"].split(", ")
        for syn in syns:
            key = SynToKey(syn)
            if key == "":
                continue
            if key in synonyms_dict and synonyms_dict[key] != compound:
                logging.warning(
                    "Warning: non-unique synonym "
                    "{}: {} or {}\n".format(syn, compound, synonyms_dict[key])
                )
                synonyms_dict[key] = compound

    return compounds_dict, synonyms_dict 


def SynToKey(syn):
    syn = re.sub(r"[^a-zA-Z0-9+-]", '', syn)
    return syn.lower()

# Second Run
# Returns two hashes:
# media => list of components
# media => attributes
# Handles both the media file and the mixes file
# Does NOT look up compound synonyms or otherwise check the results
def LoadMedia(metadir, compounds, synonyms, unknownComponents, reuseComponents):
    # Uses subfunctions ParseMediaFile, ExpandMedia, SetupComponentList,
    # SetupExclude
    mediaFile = os.path.join(metadir, "media")
    mixFile = os.path.join(metadir, "mixes")
    if not os.path.isfile(mediaFile):
        raise Exception("No such file: " + mediaFile + "\n")
    if not os.path.isfile(mixFile):
        raise Exception("No such file: " + mixFile + "\n")

    # Init mediaExclude - media => excluded compound => 1
    # This is first occurence of mediaExclude
    mediaExclude = {}
    # The following are dicts
    # media is key => list<list of len(3)>
    media, mediaAttr, mediaExclude = ParseMediaFile(mediaFile, mediaExclude)
    mix, mixAttr, mediaExclude = ParseMediaFile(mixFile, mediaExclude)

    # Validation:
    # Each media must have a Description
    # media should not also be mixes
    for k in mediaAttr.keys():
        med, attr = k, mediaAttr[k]
        if not ("Description" in attr and attr["Description"] != ""):
            raise Exception("No description for media " + med)
        if med in mix:
            raise Exception("Media {} is also a mix".format(med))

    # Each mix must have a Description and a numeric X value
    # mixes should not also be media
    for k in mixAttr.keys():
        mx, attr = k, mixAttr[k]
        if not ("Description" in attr and attr["Description"] != ""):
            raise Exception("No description for mix " + mx)
        if not ("X" in attr and (re.search(r"^[0-9]+[.]?[0-9]*$", attr["X"]))):
            raise Exception("Invalid X for mix " + mx)
        if mx in media:
            raise Exception("Mix {} is also a media".format(mx))

    # Expand media in terms of other media
    nMaxCycles = 100
    for nCycle in range(nMaxCycles + 2):
        if nCycle >= nMaxCycles:
            raise Exception(
                "Too many cycles of media expansion -- is there " "a cycle?\n"
            )
        nChanges = 0
        for md in media.keys():
            nChanges += ExpandMedia(md, media, mediaExclude)[0]
        if nChanges == 0:
            break

    # Replace compound synonyms with compounds, and record any that are not
    # known or are duplicates
    # Below media[k] is a list of length 3
    for k in media.keys():
        media, unknownComponents, reuseComponents = SetupComponentList(
            k,
            media[k],
            media,
            mix,
            unknownComponents,
            reuseComponents,
            compounds,
            synonyms,
        )
    for k in mix.keys():
        media, unknownComponents, reuseComponents = SetupComponentList(
            k,
            mix[k],
            media,
            mix,
            unknownComponents,
            reuseComponents,
            compounds,
            synonyms,
        )

    # Replace compound synonyms with compounds
    # Warn if excluded compound is in the media (and remove it)
    # Define "minus" mixes as needed
    for k in mediaExclude.keys():
        media, mix, mixAttr, unknownComponents = SetupExclude(
            k,
            mediaExclude[k],
            unknownComponents,
            media,
            mix,
            mixAttr,
            compounds,
            synonyms,
        )
    results_dict = {
            "media_dict": media,
            "mix_dict": mix,
            "mixAttr": mixAttr,
            "compounds_dict": compounds,
            "synonyms_dict": synonyms
            }
    return results_dict


def SetupExclude(
    media, excludeHash, unknownComponents, media_dict, mix, mixAttr, compounds, synonyms
):
    # Uses subfunctions FindCompound,
    excluded = []
    for orig in excludeHash.keys():
        compound = FindCompound(orig, compounds, synonyms)
        if compound is not None:
            excluded.append(compound)
        else:
            unknownComponents[orig] = 1
    # Dict of excluded compound => number of alterations
    excluded_dict = {x: 0 for x in excluded}

    # And update the media and any incorporated mixes to actually exclude
    updatedComponents = []
    for component in media_dict[media]:
        compound, number, units = component
        if compound in excluded_dict:
            excluded_dict[compound] += 1
        elif compound in mix:
            # compound is a mix -- check if it needs to be altered
            mixExcluded = {}
            mixKeep = []
            for mixComponent in mix[compound]:
                mixCompound = mixComponent[0]
                if mixCompound in excluded_dict:
                    mixExcluded[mixCompound] = 1
                else:
                    mixKeep.append(mixComponent)
            if len(mixExcluded.keys()) > 0:
                # Define the new mix and use it instead
                minusString = " ".join(
                    ["minus " + x for x in sorted(mixExcluded.keys())]
                )
                mixNew = compound + " " + minusString
                if mixNew in mix:
                    if not (len(mix[mixNew])) == (len(mixKeep)):
                        raise Exception(
                            "Wrong number of components for mix "
                            + mixNew
                            + " which can also be built via exclude"
                            " from " + compound
                        )
                else:
                    mix[mixNew] = mixKeep
                    mixAttr[mixNew] = mixAttr[compound]
                    mixAttr[mixNew]["Description"] += " " + minusString
                updatedComponents.append([mixNew, number, units])
                for mixCompound in mixExcluded.keys():
                    excluded_dict[mixCompound] += 1
            else:
                updatedComponents.append(component)
        else:
            updatedComponents.append(component)
    media_dict[media] = updatedComponents
    for k in excluded_dict.keys():
        compound, count = k, excluded_dict[k]
        if count is None:
            raise Exception("B112 " + media + compound)
        if not count > 0:
            logging.warning(
                "Warning: excluding {} from media".format(compound)
                + "{} had no effect\n".format(media)
            )

    return [media_dict, mix, mixAttr, unknownComponents]


def ParseMediaFile(mediaFile, mediaExclude):
    # Uses subfunction ListValidUnits
    comp = {} # returns as "media or mix"
    attr = {} # returns as "mediaAttr or mixAttr"
    COMPOUND, NUMBER, UNITS = 0, 1, 2
    validUnits = {x: 1 for x in ListValidUnits()}
    validAttr = {x: 1 for x in ["Description", "Minimal", "X"]}
    curMedia = None
    readingCompounds = 0
    with open(mediaFile, "r") as f:
        mf_lines = f.read().split("\n")
    for j in range(len(mf_lines)):
        l = mf_lines[j]
        # Handle DOS mode files
        l = re.sub(r'[\r\n]+$', '', l)
        # strip trailing fields that are empty
        l = re.sub(r'\t+$', '', l)
        F = l.split("\t")
        if len(F) == 1:
            curMedia = None
        elif re.search(r"^#", F[0]):
            pass
        elif len(F) == 2:
            att, val = F
            if att == "Media":
                curMedia = F[1]
                curMedia = re.sub(r" +$", "", curMedia)
                if curMedia in comp:
                    raise Exception("Duplicate media entry for " + curMedia)
                comp[curMedia] = []
                attr[curMedia] = {}
                readingCompounds = 0
            elif att in validAttr:
                if curMedia is None:
                    raise Exception(
                        "No media id yet at line:\n" "{}\n in {}".format(l, mediaFile)
                    )
                if att in attr[curMedia]:
                    raise Exception(
                        "Duplicate attr {} for media {}".format(att, curMedia)
                    )
                attr[curMedia][att] = val
            else:
                raise Exception("Invalid media attribute {}".format(F[0]))
        elif len(F) == 3:
            # LINE 324 in Perl
            if curMedia is None:
                raise Exception(
                    "No media id yet at line:\n" "{}\nin {}".format(
                        l, mediaFile)
                    )
            if (
                (re.search(r"^Controlled", F[0]))
                and F[1] == "Concentration"
                and F[2] == "Units"
            ):
                readingCompounds = 1
            else:
                if not readingCompounds == 1:
                    raise Exception(
                        "No compounds header for " "{} at \n{}\n".format(curMedia, l)
                    )
                compound, concentration, units = F
                compound = re.sub(r" +$", "", compound)
                concentration = re.sub(r"^ +", "", concentration)
                concentration = re.sub(r" +$", "", concentration)
                if concentration == "-" and units == "-":
                    if curMedia in mediaExclude:
                        mediaExclude[curMedia][compound] = 1
                    else:
                        mediaExclude[curMedia] = {compound: 1}
                else:
                    cr = concentration
                    if not (
                        cr == ""
                        or re.search(r"^\d+$", cr)
                        or (re.search(r"^\d+[.]\d*$", cr))
                        or (re.search(r"^\d+[.]?\d*[eE][+-]\d+$", cr))
                    ):
                        raise Exception(
                            "Invalid concentration "
                            + concentration
                            + " in line\n"
                            + l
                            + "\nfor "
                            + curMedia
                            + "\n"
                            "in " + mediaFile
                        )
                    units.strip()
                    if units == "":
                        raise Exception("B111: " + l)
                    if units not in validUnits:
                        raise Exception(
                            "Invalid unit "
                            + units
                            + " for\n"
                            + l
                            + "\nin "
                            + curMedia
                            + ", "
                            + mediaFile
                        )
                    comp[curMedia].append([compound, cr, units])
        else:
            raise Exception("Wrong number of fields in:\n"
                    "{}\nIn file {} at line {} ".format(l, mediaFile, j) \
                            + " Num Fields: {}".format(len(F)))
    return [comp, attr, mediaExclude]


def ExpandMedia(media_str, media_dict, mediaExclude):
    # Uses no sub functions
    components = []
    nExpand = 0
    for row in media_dict[media_str]:
        comp, conc, units = row
        conc = float(conc)
        if comp == media_str:
            raise Exception(
                "Invalid definition of medium "
                "{} -- it includes itself\n".format(media_str)
            )
        if comp in media_dict:
            # comp is another media
            if comp in mediaExclude:
                raise Exception(
                    "Invalid definition of medium "
                    "{} - it includes medium {} with an exclude".format(media_str, comp)
                    + "compound"
                )
            if not units == "X":
                raise Exception(
                    "Invalid definition of medium "
                    "{} - it includes medium {} but units are not X".format(
                        media_str, comp
                    )
                )
            for subComponent in media_dict[comp]:
                comp2, conc2, units2 = subComponent
                conc2 = float(conc2)
                conc2 *= conc
                components.append([comp2, conc2, units2])
            nExpand += 1
        else:
            components.append(row)
    media_dict[media_str] = components
    return [nExpand, mediaExclude]

# media is a str, l is a list of length 3, everything else dict? 
def SetupComponentList(
    media,
    l,
    media_dict,
    mix_dict,
    unknownComponents,
    reuseComponents,
    compounds,
    synonyms,
):
    # Uses subfunctions FindCompound
    # 1 or 0
    isMedia = True if media in media_dict else False
    isMix = True if media in mix_dict else False
    if not (isMedia or isMix):
        raise Exception("Unknown media:\n"
                "{}".format(media))
    COMPOUND = 0
    # transfer synonyms or record that it is unknown
    for row in l:
        orig, undef, units = row
        origInMix = 1 if orig in mix_dict else 0
        if isMedia == 1 and origInMix == 1:
            if not units == "X":
                raise Exception("Mix must be included with X units")
        else:
            compound = FindCompound(orig, compounds, synonyms)
            if compound is not None:
                row[COMPOUND] = compound
            else:
                unknownComponents[row[COMPOUND]] = 1
    # Record repeat entries
    seen = {}
    for row in l:
        compound = row[COMPOUND]
        if compound in seen:
            if compound in reuseComponents:
                reuseComponents[compound][media] = 1
            else:
                reuseComponents[compound] = {media:1}
        seen[compound] = 1

    return [media_dict , unknownComponents, reuseComponents]


# Complete - syn is str, compounds and synonyms are dicts 
# Run 4th - after GetMediaComponents
def FindCompound(syn, compounds, synonyms):
    # Uses subfunction SynToKey
    if syn in compounds:
        return syn
    key = SynToKey(syn)
    if key in synonyms:
        return synonyms[key]
    return None

# Run 3rd
def GetMediaComponents(media_str, media_dict, mixAttr, mix):
    """
    returns undef for unknown media; otherwise, a  
    list of [compound,concentration,units,mix]
    where mix is empty (not undefined) unless the compound was 
    included indirectly via a mix
    """
    if media_str not in media_dict:
        return None
    out = []
    for row in media_dict[media_str]:
        # comp, units are str. conc float/int
        comp, conc, units = row
        conc = float(conc)
        if comp in mix:
            if not units == "X":
                raise Exception("Units for mix {} in media {} are not X".format(
                    comp, media_str))
            if not "X" in mixAttr[comp]:
                raise Exception("No X value for mix " + comp)
            rel = conc/float(mixAttr[comp]["X"])
            if not rel > 0 and rel < 1e4:
                raise Exception("Invalid relative X {} for {} in {}".format(
                    rel, comp, media_str))
            for row2 in mix[comp]:
                comp2, conc2, units2 = row2
                conc2 = float(conc2)
                out.append([comp2, conc2*rel, units2, comp])
        else:
            out.append([comp, conc, units, ""])

    return out

# Run last (5th)
def GetMixComponents(mix, mix_dict):
    if mix not in mix_dict:
        return None
    else:
        return mix_dict[mix]




def ListValidUnits():
    return "g/L mg/L ug/L M mM uM vol% ml/L X".split(" ")
