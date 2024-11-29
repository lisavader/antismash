from collections import OrderedDict

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup, docs_link
from antismash.common.layers import OptionsLayer, RegionLayer, RecordLayer

from .results import TerpeneResults
from .terpene_analysis import _load_hmm_properties

_hmm_properties = _load_hmm_properties()

def will_handle(_products: list[str], categories: set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return "terpene" in categories

def get_domain_description(prediction) -> str:
    description = _hmm_properties[prediction.type].description
    return description

def format_subtype(prediction) -> str:
    if not prediction.subtypes:
        if _hmm_properties[prediction.type].subtypes:
            return "unknown"
        return "none"
    descriptions = []
    for subtype in prediction.subtypes:
        try:
            short_description = "".join(_hmm_properties[subtype].description.split("; ")[1:])
        except:
            short_description = subtype.description
        descriptions.append(short_description)
    return " or ".join(descriptions)

def format_reactions(prediction) -> Markup:
    markup = Markup()
    for reaction in prediction.reactions:
        substrates = [substrate.name for substrate in reaction.substrates]
        products = [product.name for product in reaction.products]
        for product in products:
            markup += Markup(f"<dd>{', '.join(substrates)} &rarr; {product}</dd>")
    return markup

def format_glossary_rows(protoclusters) -> Markup:
    markup = Markup()
    name_mappings = {}
    for protocluster in protoclusters:
        for domains in protocluster["prediction"].cds_predictions.values():
            for domain in domains:
                for reaction in domain.reactions:
                    compounds = reaction.substrates + reaction.products
                    for compound in compounds:
                        if compound.extended_name:
                            name_mappings[compound.name] = compound.extended_name
    ordered_name_mappings = OrderedDict(sorted(name_mappings.items()))
    for name, extended_name in ordered_name_mappings.items():
        markup += Markup("<tr>"
                         f"<td>{name}</td>"
                         f"<td>{extended_name.capitalize()}</td>"
                         "</tr>")
    return markup

def format_domain_types(domain_preds) -> str:
    types = []
    for domain_pred in domain_preds:
        types.append(domain_pred.type)
    return " + ".join(types)

def generate_html(region_layer: RegionLayer, results: TerpeneResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML output for the module """
    html = HTMLSections("terpene")

    protoclusters = []
    for cluster in region_layer.get_unique_protoclusters():
        if cluster.product_category == "terpene":
            protoclusters.append({
                "cluster" : cluster,
                "prediction" : results.cluster_predictions[cluster.get_protocluster_number()]})

    details_template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    sidepanel_template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))

    tooltip = "Some tooltip information for the details panel."
    details = details_template.render(protoclusters = protoclusters, tooltip = tooltip,
                                      get_domain_description = get_domain_description,
                                      format_subtype = format_subtype, format_reactions = format_reactions,
                                      format_domain_types = format_domain_types)
    html.add_detail_section("Terpene", details, class_name="terpene")

    tooltip = "Some tooltip information for the sidepanel."
    sidepanel = sidepanel_template.render(protoclusters = protoclusters, tooltip = tooltip,
                                          format_glossary_rows = format_glossary_rows)
    html.add_sidepanel_section("Terpene", sidepanel, class_name="terpene")

    return html