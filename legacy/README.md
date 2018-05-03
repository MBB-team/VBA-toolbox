# Legacy code

This folder contains deprecated functions. To ensure backward compatibility, calling any of those functions will try to fallback to a replacement function, or produce an error if no alternative can be found.

In any case, those functions also issue some warning to help adapting user code.
