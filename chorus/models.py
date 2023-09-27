import uuid

from django.contrib.auth.models import User
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
        return f"{self.protein.name} {self.original}{self.position}{self.mutated}"

    class Meta:
        ordering = ["id"]

class ChorusSession(models.Model):
    """
    class for storing user submitted session data that contains filter condition
    """
    user = models.ForeignKey(User, on_delete=models.RESTRICT, related_name="chorus_session", null=True, blank=True)
    link_id = models.UUIDField(default=uuid.uuid4, editable=False)
    description = models.TextField()
    file = models.FileField(upload_to='chorus/files/', blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.link_id

    class Meta:
        ordering = ["id"]