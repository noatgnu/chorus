from django.db import models

class Protein(models.Model):
    """class for storing protein name to be used later as key for variant data class"""
    name = models.CharField(max_length=100)
    description = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ["id"]

class Variant(models.Model):
    """class for storing variant data that use protein as key and has 3 other columns storing position, original residue and mutated residue"""
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    position = models.IntegerField()
    original = models.CharField(max_length=100)
    mutated = models.CharField(max_length=100)
    score = models.FloatField()
    pathogenicity = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.protein

    class Meta:
        ordering = ["id"]

