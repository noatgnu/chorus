

class VariantProteinRouter:
    """
    A router to control all database operations for the model Variant.
    """
    route_app_labels = {'variant_data', 'protein_data'}
    def db_for_read(self, model, **hints):
        """
        Attempts to read protein and variant models go to protein_db.
        """
        if model._meta.app_label in self.route_app_labels:
            return 'variant'
        return None

    def db_for_write(self, model, **hints):
        """
        Attempts to write protein and variant models go to protein_db.
        """
        if model._meta.app_label in self.route_app_labels:
            return 'variant'
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the variant_data app is involved.
        """
        if (obj1._meta.app_label in self.route_app_labels) or \
           (obj2._meta.app_label in self.route_app_labels):
           return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Make sure the variant_data app only appears in the 'variant'
        database.
        """
        if app_label in self.route_app_labels:
            return db == 'variant'
        return None
