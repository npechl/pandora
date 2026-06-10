from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field


Severity = Literal["warning", "error"]
ResultStatus = Literal["success", "warning", "failed"]


class Diagnostic(BaseModel):
    code: str
    severity: Severity
    message: str
    entry_id: str | None = None
    field: str | None = None
    context: dict[str, Any] | None = None


class DiagnosticBundle(BaseModel):
    warnings: list[Diagnostic] = Field(default_factory=list)
    errors: list[Diagnostic] = Field(default_factory=list)


class AppliedPolicyRef(BaseModel):
    """Compact reference to a policy that was applied at a pipeline stage."""
    policy_id: str
    policy_name: str
    policy_version: str
