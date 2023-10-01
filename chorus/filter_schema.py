import six
from filters.schema import base_query_params_schema
from filters.validations import IntegerLike

protein_query_schema = base_query_params_schema.extend(
    {
        "id": IntegerLike(),
        "name": six.text_type,
    }
)

variant_query_schema = base_query_params_schema.extend(
    {
        "id": IntegerLike(),
        "protein": six.text_type,
        "position": IntegerLike(),
        "original": six.text_type,
        "mutated": six.text_type,
        "pathogenicity": six.text_type,
    }
)

chorus_session_query_schema = base_query_params_schema.extend(
    {
        "id": IntegerLike(),
        "user": IntegerLike(),
        "link_id": six.text_type,
    }
)