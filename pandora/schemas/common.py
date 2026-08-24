from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

Severity = Literal["warning", "error"]
ResultStatus = Literal["success", "warning", "failed"]


class Diagnostic(BaseModel):
    """One warning or error raised during parsing/canonicalisation, with
    its code, message, and context.

    Attributes:
        code: A short machine-readable code identifying the issue
            (e.g. "RESIDUE_NUMBER_COLLISION").
        severity: "warning" or "error".
        message: A human-readable description of the issue.
        entry_id: The entry this diagnostic applies to, if known.
        field: The field/category this diagnostic relates to, if
            applicable.
        context: Additional structured details about the issue.
    """

    code: str
    severity: Severity
    message: str
    entry_id: str | None = None
    field: str | None = None
    context: dict[str, Any] | None = None


class DiagnosticBundle(BaseModel):
    """Every warning and error collected during one parsing/
    canonicalisation run.

    Attributes:
        warnings: Every warning-level diagnostic collected.
        errors: Every error-level diagnostic collected.
    """

    warnings: list[Diagnostic] = Field(default_factory=list)
    errors: list[Diagnostic] = Field(default_factory=list)


class AppliedPolicyRef(BaseModel):
    """Compact reference to a policy that was applied at a pipeline stage.

    Attributes:
        policy_id: The policy's unique identifier.
        policy_name: The policy's human-readable name.
        policy_version: The policy's version string.
    """

    policy_id: str
    policy_name: str
    policy_version: str
