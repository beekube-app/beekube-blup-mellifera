"""Catalogue des caractères standard (spec volet 3 §2), figé dans le contrat et donc dans le SDK.

Un critère Beekube déclare son caractère (`CriterionIn.trait`) ; le moteur s'en sert pour rattacher les notes
Index Mellifera des ancêtres résolus (colonne `mellifera`, noms de `mellifera/traits.csv`). Les colonnes
BeeBreed et Karl Kehrle sont documentaires, pour les exports futurs de Beekube. Correspondances une pour une.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional


@dataclass(frozen=True)
class Trait:
    label: str
    kind: Literal["note", "percentage", "harvest"]
    mellifera: Optional[str] = None      # trait.name de la base du crawler (mellifera/traits.csv)
    beebreed: Optional[str] = None       # caractère publié par BeeBreed
    kehrle: Optional[str] = None         # champ du registre Karl Kehrle
    direction: Literal[1, -1] = 1        # 1 : valeur haute = mieux


STANDARD_TRAITS: dict[str, Trait] = {
    "honey_yield": Trait("Récolte de miel (kg)", "harvest", "honeyAmt", "Honey yield", "honeykgpromean"),
    "honey_yield_average": Trait("Récolte de miel moyenne (kg)", "harvest", "aveHonAmt"),
    "honey_score": Trait("Indice miel (note)", "note", "honeyIdx", None, "honeyindex"),
    "average_crop_score": Trait("Indice de récolte moyen (note)", "note", "aveCropIdx"),
    "spring_crop": Trait("Récolte de printemps (note)", "note", "sprCrop", None, "cropearly"),
    "mid_season_crop": Trait("Récolte de mi-saison (note)", "note", "midCrop"),
    "summer_crop": Trait("Récolte d'été (note)", "note", "sumCrop", None, "cropsumm"),
    "autumn_crop": Trait("Récolte d'automne (note)", "note", "falCrop", None, "croplate"),
    "pollen_crop": Trait("Récolte de pollen (note)", "note", "polCrop"),
    "gentleness": Trait("Douceur", "note", "mildness", "Defensive behavior", "temperofbees"),
    "calmness": Trait("Tenue au cadre", "note", "framHold", "Steadiness on comb", "calmbehaviour"),
    "swarming": Trait("Essaimage (note haute = essaime peu)", "note", "swarming", "Swarming drive", "swarmtendency"),
    "overwintering": Trait("Hivernage", "note", "winter", "Robustness in winter", "overwintering"),
    "spring_development": Trait("Développement de printemps", "note", None, "Development in Spring"),
    "colony_strength": Trait("Force de la colonie", "note", None, "Colony strength"),
    "brood_vitality": Trait("Vitalité du couvain", "note", "broVital", None, "vitalityofbrood"),
    "bee_vitality": Trait("Vitalité des abeilles", "note", "beeVital", None, "vitalityofbees"),
    "fertility": Trait("Fécondité", "note", "fertil", None, "fertility"),
    "comb_building": Trait("Bâtisse", "note", "bldCombs", None, "buildingcombs"),
    "propolis": Trait("Propolis", "note", "propolis", None, "buildingpropol"),
    "hygienic_24h": Trait("Test d'aiguille, couvain nettoyé à 24 h (%)", "percentage", "hyg24h",
                          "rate of opened cells", "hygprozafter24h"),
    "hygienic_48h": Trait("Test d'aiguille, couvain nettoyé à 48 h (%)", "percentage", "hyg48h", None,
                          "hygienicprozafter48h"),
    "hygienic_score": Trait("Indice de comportement hygiénique (note)", "note", "hygBehaIdx", None, "hygienicindex"),
    "varroa_score": Trait("Indice varroa (note)", "note", "varIdx", "Varroa index"),
    "varroa_infestation": Trait("Infestation varroa (mesure)", "percentage", None,
                                "Varroa infestation development", "varroanumber", direction=-1),
    "smr": Trait("SMR / VSH (%)", "percentage", None, "SMR", "varroavsh"),
    "recapping": Trait("Recapping (%)", "percentage", None, "Recapping"),
    "chalkbrood": Trait("Couvain plâtré (note haute = résistant)", "note", None, "Chalkbrood"),
    "nosema": Trait("Nosémose (note haute = résistant)", "note", None, "Nosemosis", "nosemaindex"),
    "cbpv": Trait("Paralysie chronique (note haute = résistant)", "note", None, "Chronic Bee Paralysis"),
}

TraitName = Literal[tuple(STANDARD_TRAITS)]


def trait_description() -> str:
    """Description OpenAPI de CriterionIn.trait : le catalogue et ses correspondances."""
    lines = ["Caractère standard du critère ; sert au moteur à rattacher les notes Index Mellifera des ancêtres "
             "résolus, et à Beekube pour ses exports. Catalogue (identifiant : libellé [type] ; Index Mellifera, "
             "BeeBreed, Karl Kehrle) :"]
    for key, t in STANDARD_TRAITS.items():
        refs = " ; ".join(f"{src} {name}" for src, name in (("Index Mellifera", t.mellifera), ("BeeBreed", t.beebreed),
                                                           ("Karl Kehrle", t.kehrle)) if name)
        sens = ", valeur haute = moins bien" if t.direction == -1 else ""
        lines.append(f"- {key} : {t.label} [{t.kind}{sens}]" + (f" ; {refs}" if refs else ""))
    return "\n".join(lines)
