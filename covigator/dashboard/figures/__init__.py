from covigator.database.model import VariantType

VARIANT_TYPE_COLOR_MAP = {
    VariantType.SNV.name: "#8da0cb",
    VariantType.INSERTION.name: "#fc8d62",
    VariantType.DELETION.name: "#66c2a5",
    VariantType.MNV.name: "#e78ac3",
    VariantType.COMPLEX.name: "#a6d854"
}
